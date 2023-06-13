/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adaptiveImplicitAdvection.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"
#include "fvc.H"
#include "fvcFluxLimit.H"
#include "CourantNoFunc.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * * * //

template<class Type>
adaptiveImplicitAdvection<Type>::adaptiveImplicitAdvection
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    convectionScheme<Type>(mesh, faceFlux),
    dict_(is),
    tCorrectionScheme_
    (
        surfaceInterpolationScheme<Type>::New
        (
            mesh, faceFlux, dict_.lookup("correctionScheme")
        )
    ),
    nCorr_(readLabel(dict_.lookup("nCorr"))),
    offCentreName_(dict_.lookup("offCentre")),
    CoLimit_(readScalar(dict_.lookup("CoLimit"))),
    FCTlimit_(dict_.lookup("FCTlimit")),
    FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
    FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0)))
{}

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

template<class Type>
const surfaceScalarField& 
adaptiveImplicitAdvection<Type>::offCentre() const
{
    return this->mesh().objectRegistry::template
         lookupObject<surfaceScalarField>(offCentreName_);
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
adaptiveImplicitAdvection<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Un flux-corrected interpolation
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            tCorrectionScheme_().interpolate(vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
adaptiveImplicitAdvection<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
adaptiveImplicitAdvection<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Matrix to be returned will only have a source term
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    tfvm.ref() += fvcDiv(faceFlux, vf);
    
    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
adaptiveImplicitAdvection<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Schemes needed
    localMax<scalar> maxInterp(mesh);
    EulerDdtScheme<Type> backwardEuler(mesh);
    gaussConvectionScheme<Type> upwindConvect
    (
        mesh,
        faceFlux,
        tmp<surfaceInterpolationScheme<Type>>(new upwind<Type>(mesh, faceFlux))
    );
    
    // Update the implicit/explicit split
    surfaceScalarField Cof(maxInterp.interpolate(CourantNo(faceFlux, dt)));
    surfaceScalarField ImEx("ImEx", 0.5*(sign(Cof - CoLimit_) + 1));
    surfaceScalarField oImEx("oImEx", offCentre()*ImEx);

    // Initialise the divergence to be the old time, upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect.fvcDiv((1-oImEx)*faceFlux, vf.oldTime())
    );
    tConvection.ref().rename("convection(" + faceFlux.name() + ',' + vf.name() + ')');

    // Create temporary field to advect and apply explicit upwind
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf.oldTime() - dt*tConvection(),
        vf.boundaryField().types()
    );

    // Implicit upwind
    if (max(oImEx).value() > SMALL)
    {
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(T)
          + fvMatrix<Type>(upwindConvect.fvmDiv(oImEx*faceFlux, T))
        );
        TEqn.solve();
        tConvection.ref() += fvc::div(TEqn.flux());
        T.oldTime() = vf.oldTime() - dt*tConvection();
    }
    else T.oldTime();

    // Correction at the old time
    surfaceScalarField fluxCorrOld = (1-offCentre())*faceFlux.oldTime()
                               *tCorrectionScheme_->correction(vf.oldTime());
    // Correction to be updated
    surfaceScalarField fluxCorr = fluxCorrOld;
    
    // HO corrections
    for(int iCorr = 0; iCorr < nCorr_; iCorr++)
    {
        fluxCorr = fluxCorrOld
                 + offCentre()*faceFlux*tCorrectionScheme_->correction(T);

        // Add HO fluxes
        if (iCorr < nCorr_-1)
        {
            T = T.oldTime() - dt*fvc::div(fluxCorr);
        }
    }
    
    // Limit the fluxes with Zalesak FCT limiter if needed
    if (FCTlimit_ && nCorr_ > 0)
    {
        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_);
        
        if (FCTmin_ < FCTmax_)
        {
            fvc::fluxLimit(fluxCorr, T, FCTmin, FCTmax, dt);
        }
        else if (max(oImEx).value() < SMALL)
        {
            fvc::fluxLimit(fluxCorr, T.oldTime(), vf.oldTime(), dt);
        }
        else
        {
            fvc::fluxLimit(fluxCorr, T.oldTime(), dt);
        }
    }
    
    if (nCorr_ > 0)
    {
        T = T.oldTime() - dt*fvc::div(fluxCorr);
        tConvection.ref() += fvc::div(fluxCorr);
    }

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
