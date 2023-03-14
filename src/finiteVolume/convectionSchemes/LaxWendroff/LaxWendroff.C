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
#include "uncorrectedSnGrad.H"
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void LaxWendroff<Type>::calculateAntiD
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& offCentre
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();
    localMax<scalar> maxInterp(this->mesh());

    // Calculate necessary additional fields for the correction

    // Gradient of T in the cell centre to cell centre direction
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    surfaceScalarField snGradT = snGrad.snGrad(vf);

    // The correction in space
    antiD() = 0.5*mag(faceFlux)*snGradT/rdelta;

    // Apply a correction in time if needed
    if (timeCorrector_ == "advective")
    {
        // The full velocity field from the flux and correct
        surfaceVectorField Uf = fvc::interpolate
        (
            fvc::reconstruct(faceFlux), "LaxWendroff_velocity"
        );
        Uf += (faceFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        // Full face gradient of T
        fv::leastSquaresGrad<Type> grad(mesh);
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        > gradT = fvc::interpolate(grad.grad(vf), "LaxWendroff_gradient");
        gradT += (snGradT - (gradT & mesh.delta())*rdelta)
             * mesh.delta()*rdelta;

        // add advetive form time correction
        antiD() -= max(1-2*offCentre, scalar(0))*0.5*faceFlux*dt*(Uf & gradT);
    }
    else if(timeCorrector_ == "flux")
    {
        antiD() -= max(1-2*offCentre, scalar(0))*0.5*faceFlux*dt
         *fvc::interpolate(fvc::div(faceFlux, vf, "LaxWendroff_div"), "LaxWendroff_idiv");
    }

    // Smooth where offCentre>0
    surfaceVectorField V("antiDV", linearInterpolate(fvc::reconstruct(antiD())));
    surfaceScalarField imp = offCentre/(offCentre+SMALL);
    imp = maxInterp.interpolate(fvc::localMax(imp));
    //imp = linearInterpolate(fvc::localMax(imp));
    antiD() = imp*(V & mesh.Sf()) + (1-imp)*antiD();

    /*// Limit to obey Courant number restriction
    volScalarField CoV("CoV", 4*CourantNo(antiD(), dt));
    surfaceScalarField Cof("Cof", maxInterp.interpolate(CoV));
    antiD() /= max(scalar(1), Cof);
    */
    volScalarField CoV("CoV", CourantNo(antiD(), dt));
    Info << "Anti-diffusive Courant number max: " << max(CoV).value()
         << endl;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
LaxWendroff<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Just use the low order interpolate for now
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

    // Calculate the local off centering for each face
    localMax<scalar> maxInterp(this->mesh());
    volScalarField C = CourantNo(faceFlux, dt);

    surfaceScalarField offCentre
    (
        IOobject("offCentre", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    );
    if (offCentre_ < 0)
    {
        volScalarField offCentreC
        (
            "offCentreC",
            max(1-1/(C+0.25), scalar(0))
        );
        offCentre = maxInterp.interpolate(offCentreC);
    }
    else {offCentre == offCentre_;}
    
    Info << "offCentre goes from " << min(offCentre).value() << " to " 
         << max(offCentre).value() << endl;
    
    // Initialise the divergence to be the first-order upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv((1-offCentre)*faceFlux, vf)
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
    if (mag(offCentre_) > SMALL)
    {
        EulerDdtScheme<Type> backwardEuler(this->mesh());
        fvMatrix<Type> fvmT
        (
            backwardEuler.fvmDdt(T)
          + upwindConvect().fvmDiv(offCentre*faceFlux, T)
        );
        fvmT.solve();
    
        // Add the low order implicit divergence
        tConvection.ref() += upwindConvect().fvcDiv(offCentre*faceFlux, T);
    }
    
    // Calculate, apply (and update) the correction
    if (nCorr_ > 0) calculateAntiD(faceFlux, T, offCentre);
    for(label iCorr = 1; iCorr < nCorr_; iCorr++)
    {
        calculateAntiD
        (
            faceFlux,
            T - dt*fvc::div(antiD()),
            offCentre
        );
    }
    if (nCorr_ > 0)
    {
        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_);
        
        // Limit the fluxes with Zalesak FCT limiter if needed
        if (FCTlimit_)
        {
            if (FCTmin_ < FCTmax_)
            {
                fvc::fluxLimit(antiD(), T, FCTmin, FCTmax, dt);
            }
            else if (mag(offCentre_) < SMALL)
            {
                fvc::fluxLimit(antiD(), T, vf, dt);
            }
            else
            {
                fvc::fluxLimit(antiD(), T, dt);
            }
        }
    
        tConvection.ref() += fvc::div(antiD());
    }
    
    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
