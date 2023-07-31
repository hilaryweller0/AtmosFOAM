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

template<class Type>
tmp<surfaceScalarField>
adaptiveImplicitAdvection<Type>::calcOffCentre
(
    const surfaceScalarField& volFlux
) const
{
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    localMax<scalar> maxInterp(mesh);

    tmp<surfaceScalarField> toff
    (
        new surfaceScalarField
        (
            IOobject("offCentre", mesh.time().timeName(), mesh),
            mesh,
            dimensionedScalar("", dimless, offCentre_)
        )
    );
    surfaceScalarField& off = toff.ref();
    
    if (offCentre_ < -SMALL)
    {
        // First set offCentre_ to Courant number (on the face)
        if (mesh.objectRegistry::template foundObject<volScalarField>("Co"))
        {
            off = maxInterp.interpolate
            (
                mesh.objectRegistry::template
                     lookupObject<volScalarField>("Co")
            );
        }
        else
        {
            off = maxInterp.interpolate(CourantNo(volFlux, dt));
        }
        // Next calculate off centering from Courant number
        //off = max(scalar(0.5), 1 - 1/(off + 0.25));
        off = max(scalar(0.5), 1 - 1/off);
    }
    
    return toff;
}

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
    nCorr_(readLabel(dict_.lookup("nCorr"))),
    offCentre_(readScalar(dict_.lookup("offCentre"))),
    CoLimit_(readScalar(dict_.lookup("CoLimit"))),
    //fullSolver_(dict_.lookup("fullSolver")),
    FCTlimit_(dict_.lookup("FCTlimit")),
    FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
    FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0)))
{}

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //


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
            corrScheme(faceFlux)->interpolate(vf)
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
    gaussConvectionScheme<Type> upwindConvectOld
    (
        mesh,
        faceFlux.oldTime(),
        tmp<surfaceInterpolationScheme<Type>>
        (
            new upwind<Type>(mesh, faceFlux.oldTime())
        )
    );
    
    // Update the implicit/explicit split
    surfaceScalarField offCentre = calcOffCentre(faceFlux);
    surfaceScalarField Cof(maxInterp.interpolate(CourantNo(faceFlux, dt)));
    surfaceScalarField ImEx("ImEx", 0.5*(sign(Cof - CoLimit_) + 1));
    surfaceScalarField oImEx("oImEx", offCentre*ImEx);
    const Switch anyImplicit = max(oImEx).value() > SMALL;

    // Create temporary field to advect
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf,
        vf.boundaryField().types()
    );
    T.oldTime() == vf.oldTime();
    
    // Accumulate the total flux, starting from the old time step terms
    surfaceScalarField totalFlux = (1-offCentre)*faceFlux.oldTime()*
    (
        upwindConvectOld.interpolate(faceFlux.oldTime(), T.oldTime())
      + corrScheme(faceFlux.oldTime())->correction(T.oldTime())
    );

    // Low and HO corrections
    for(int iCorr = 0; iCorr < nCorr_; iCorr++)
    {
        surfaceScalarField totalFluxNew = totalFlux
            + offCentre*faceFlux*
            (
                (1-ImEx)*upwindConvect.interpolate(faceFlux, T)
              + corrScheme(faceFlux)->correction(T)
            );
    
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(T)
          + fvc::div(totalFluxNew)
          + upwindConvect.fvmDiv(oImEx*faceFlux, T)
        );
        TEqn.solve();
        if (iCorr == nCorr_-1) totalFlux = totalFluxNew + TEqn.flux();
    }

    // Limit the fluxes with Zalesak FCT limiter
    if (FCTlimit_)
    {
        T == T.oldTime();
        
        // Calculate the low order solution
        surfaceScalarField lowFlux = (1-offCentre)*faceFlux.oldTime()*
            upwindConvectOld.interpolate(faceFlux.oldTime(), T.oldTime())
          + offCentre*faceFlux*(1-ImEx)*upwindConvect.interpolate(faceFlux, T);
        
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(T)
          + fvc::div(lowFlux)
          + upwindConvect.fvmDiv(oImEx*faceFlux, T)
        );
        TEqn.solve();
        lowFlux += TEqn.flux();
        //T = T.oldTime() - dt*fvc::div(lowFlux);
        /*Info << "Low order FCT solutions goes from "
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        volScalarField Tcons = T.oldTime() - dt*fvc::div(lowFlux);
        Info << "Conservative low order FCT solutions goes from "
             << min(Tcons.internalField()).value() << " to "
             << max(Tcons.internalField()).value() << endl;*/
        
        surfaceScalarField fluxCorr = totalFlux - lowFlux;

        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_);
        
        if (FCTmin_ < FCTmax_)
        {
            fvc::fluxLimit(fluxCorr, T, FCTmin, FCTmax, dt);
        }
        else if (!anyImplicit)
        {
            fvc::fluxLimit(fluxCorr, T, vf.oldTime(), dt);
        }
        else
        {
            fvc::fluxLimit(fluxCorr, T, dt);
        }
        totalFlux = lowFlux + fluxCorr;
    }

    // Calculate the implied divergence and return
    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            "div("+faceFlux.name()+','+vf.name()+')',
            fvc::div(totalFlux)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
