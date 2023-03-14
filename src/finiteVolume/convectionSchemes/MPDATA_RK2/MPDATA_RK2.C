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

#include "MPDATA_RK2.H"
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
void MPDATA_RK2<Type>::calculateAnteD
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Calculate necessary additional fields for the correction

    // The volume field interpolated onto faces for the denominator
    GeometricField<Type, fvsPatchField, surfaceMesh> Tf
         = fvc::interpolate(vf, "MPDATA_denom");

    // Stabilisation
    Tf += dimensionedScalar("", Tf.dimensions(), gauge_ + SMALL);

    // Gradient of T in the cell centre to cell centre direction
    fv::uncorrectedSnGrad<Type> snGrad(mesh);
    surfaceScalarField snGradT = snGrad.snGrad(vf);

    // The correction in space
    anteD() = 0.5/Tf*mag(faceFlux)*snGradT/rdelta;
    
    // Limit to obey Courant number restriction
    anteD() = sign(anteD())*min(mag(anteD()), faceVol_/(4*dt));

/*    // Write out the ante-diffusive velocity if needed
    if (mesh.time().writeTime())
    {
        // Ante-diffusive velocity and divergence
        surfaceVectorField V
        (
            "anteDV",
            linearInterpolate(fvc::reconstruct(anteD()))
        );
        V += (anteD() - (V & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        V.write();
        
        // divergence of V
        volScalarField divAnteD("divAnteD", fvc::div(anteD()));
        divAnteD.write();
    }*/
    volScalarField CoV("CoV", CourantNo(anteD(), dt));
    Info << "Ante-diffusive Courant number max: " << max(CoV).value() << endl;

}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
MPDATA_RK2<Type>::interpolate
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
MPDATA_RK2<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
MPDATA_RK2<Type>::fvmDiv
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
MPDATA_RK2<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();

    // Initialise the divergence to be the first-order upwind divergence
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        upwindConvect().fvcDiv(faceFlux, vf)
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );
    
    // Create temporary field to advect with the temporary divergence field
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf - dt*tConvection(),
        vf.boundaryField().types()
    );
    
    calculateAnteD(faceFlux, T);
    T -= dt*anteDConvect().fvcDiv(anteD(), T);
    
    // Update the RK2 divergence
    tConvection.ref() = anteDConvect().fvcDiv(anteD(), T) + 0.5*
    (
        tConvection()
      + upwindConvect().fvcDiv(faceFlux, T)
    );
    
    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
