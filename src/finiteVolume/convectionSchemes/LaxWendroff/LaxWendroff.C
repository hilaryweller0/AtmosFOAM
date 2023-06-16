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

#include "LaxWendroff.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"
#include "fvc.H"
#include "orthogonalSnGrad.H"
#include "CourantNoFunc.H"
#include "localMax.H"
#include "leastSquaresGrad.H"
#include "fvcFluxLimit.H"
#include "fvcLocalMinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
 LaxWendroff<Type>::LaxWendroff
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
    timeCorrector_(dict_.lookupOrDefault<word>("timeCorrector", "advective")),
    FCTlimit_(dict_.lookup("FCTlimit")),
    FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
    FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0))),
    tupwindConvection_
    (
        tmp<gaussConvectionScheme<Type>>
        (
            new gaussConvectionScheme<Type>
            (
                mesh,
                faceFlux,
                tmp<surfaceInterpolationScheme<Type>>
                (
                    new upwind<Type>(mesh, faceFlux)
                )
            )
        )
    )
{
    if (offCentre_ >= 0.5 - SMALL) timeCorrector_ = "none";

    else if
    (
        timeCorrector_ != "advective"
     && timeCorrector_ != "flux"
     && timeCorrector_ != "none"
    )
    {
        FatalErrorIn("LaxWendroff::LaxWendroff")
            << "timeCorrector must be one of advective, flux or none but "
            << timeCorrector_ << " was given" << exit(FatalError);
        
    }
}

// * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

template<class Type>
void LaxWendroff<Type>::calculateFluxCorr
(
    surfaceScalarField& fluxCorr,
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& offCentre
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Schemes needed
    fv::orthogonalSnGrad<Type> snGrad(mesh);
    fv::gaussGrad<Type> grad(mesh);
    localMax<scalar> maxInterp(mesh);
    upwind<vector> upwindV(mesh, faceFlux);

    // Gradient of T in the cell centre to cell centre direction
    surfaceScalarField snGradT = snGrad.snGrad(vf);

    // Fields needed
    // Upwind cell centre (for each face)
    surfaceVectorField Cup(upwindV.interpolate(mesh.C()));
    // Full face gradient of T: average between face and upwind location
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    > gradT = linearInterpolate(grad.grad(vf));
    gradT += (snGradT - (gradT & mesh.delta())*rdelta)
         * mesh.delta()*rdelta;
    gradT = 0.5*(gradT + upwindV.interpolate(grad.grad(vf)));
    
    // The correction in space
    //fluxCorr = 0.5*mag(faceFlux)*snGradT/rdelta;
    fluxCorr = faceFlux*(gradT & (mesh.Cf() - Cup));

    // Apply a correction in time if needed
    if (timeCorrector_ == "advective")
    {
        // The full velocity field from the flux and correct
        surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(faceFlux));
        Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        // add advetive form time correction
        fluxCorr -= max(1-2*offCentre, scalar(0))*0.5*faceFlux*dt*(Uf & gradT);
    }
    else if(timeCorrector_ == "flux")
    {
        fluxCorr -= max(1-2*offCentre, scalar(0))*0.5*faceFlux*dt
         *linearInterpolate(fvc::div(faceFlux*linearInterpolate(vf)));
    }

    // Reduce the correction where offCentre > 0
    //fluxCorr *= (1 - offCentre);

    // Smooth where offCentre>0
    surfaceVectorField V(linearInterpolate(fvc::reconstruct(fluxCorr)));
    surfaceScalarField imp = offCentre/(offCentre+SMALL);
    imp = maxInterp.interpolate(fvc::localMax(imp));
    fluxCorr = imp*(V & mesh.Sf()) + (1-imp)*fluxCorr;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
LaxWendroff<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Just use the low order interpolate
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            upwindConvect().interpolate(faceFlux, vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
LaxWendroff<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
LaxWendroff<Type>::fvmDiv
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
LaxWendroff<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();

    // Schemes needed
    EulerDdtScheme<Type> backwardEuler(this->mesh());
    localMax<scalar> maxInterp(this->mesh());

    // Calculate the local off centering for each face
    surfaceScalarField offCentre
    (
        IOobject("offCentre", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, offCentre_)
    );
    if (offCentre_ < 0)
    {
        offCentre = maxInterp.interpolate
        (
            max(scalar(0), 1 - 1/(CourantNo(faceFlux, dt)+ 0.25))
        );
    }
    const Switch anyImplicit = max(offCentre).value() > SMALL;
    
    // Initialise the divergence to be the first-order upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv((1-offCentre)*faceFlux, vf.oldTime())
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );
    
    // Create temporary field to advect and the temporary divergence field
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf - dt*tConvection(),
        vf.boundaryField().types()
    );
    
    // Calculte the implicit part of the advection
    if (anyImplicit)
    {
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(T)
          + upwindConvect().fvmDiv(offCentre*faceFlux, T)
        );
        TEqn.solve();
    
        // Add the low order implicit divergence
        tConvection.ref() += fvc::div(TEqn.flux());
    }
    
    // Calculate, apply (and update) the correction
    surfaceScalarField fluxCorr
    (
        IOobject("fluxCorr", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(T.dimensions()*faceFlux.dimensions(), scalar(0))
    );
    if (nCorr_ > 0) calculateFluxCorr(fluxCorr, faceFlux, T, offCentre);
    for(label iCorr = 1; iCorr < nCorr_; iCorr++)
    {
        calculateFluxCorr
        (
            fluxCorr,
            faceFlux,
            T - dt*fvc::div(fluxCorr),
            offCentre
        );
    }
     // Limit the fluxes with Zalesak FCT limiter
    if (FCTlimit_ && nCorr_ > 0)
    {
        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_);
        
        if (FCTmin_ < FCTmax_)
        {
            fvc::fluxLimit(fluxCorr, T, FCTmin, FCTmax, dt);
        }
        else if (!anyImplicit)
        {
            fvc::fluxLimit(fluxCorr, T, vf, dt);
        }
        else
        {
            fvc::fluxLimit(fluxCorr, T, dt);
        }
    }
    
    if (nCorr_ > 0)
    {
        tConvection.ref() += fvc::div(fluxCorr);
    }

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
