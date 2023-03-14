/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "faceToPointReconstruct.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, pointPatchField, pointMesh
    >
> faceToPointReconstruct
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, pointPatchField, pointMesh> > treconField
    (
        new GeometricField<GradType, pointPatchField, pointMesh>
        (
            IOobject("volIntegrate("+ssf.name()+')', ssf.instance(), mesh),
            pointMesh(mesh),
            dimensioned<GradType>("", ssf.dimensions(), GradType::zero),
            "zeroGradient"
        )
    );
    GeometricField<GradType, pointPatchField, pointMesh>& grad = treconField();
    pointScalarField pointVol
    (
        IOobject("pointVol", ssf.instance(), mesh),
        pointMesh(mesh),
        dimensionedScalar("", dimVol, scalar(0))
    );
    
    // Loop over all faces and add contribution of gradient to the points
    forAll(ssf, faceI)
    {
        // Find all points of the face and loop over them
        const labelList& face = mesh.faces()[faceI];
        forAll(face, ip)
        {
            label iPt = face[ip];
            grad[iPt] += ssf[faceI] * mesh.Sf()[faceI]
                         / mesh.deltaCoeffs()[faceI];
            pointVol[iPt] += mesh.magSf()[faceI]/mesh.deltaCoeffs()[faceI];
        }
    }
    // And all boundary faces
    forAll(ssf.boundaryField(), patchI)
    {
        const scalarField& ssfb = ssf.boundaryField()[patchI];
        forAll(ssfb, faceI)
        {
            const labelList& face = mesh.boundaryMesh()[patchI][faceI];
            forAll(face, ip)
            {
                label iPt = face[ip];
                grad[iPt] += ssfb[faceI]
                           * mesh.boundary()[patchI].Sf()[faceI]
                          *.5/ mesh.boundary()[patchI].deltaCoeffs()[faceI];
                pointVol[iPt] += mesh.boundary()[patchI].magSf()[faceI]
                             *0.5/ mesh.boundary()[patchI].deltaCoeffs()[faceI];
            }
        }
    }
    
    grad.internalField() /= 0.5*pointVol.internalField();
    
    return treconField;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, pointPatchField, pointMesh
    >
>
faceToPointReconstruct
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, pointPatchField, pointMesh> >
         tf(fvc::faceToPointReconstruct(tssf()));
    tssf.clear();
    return tf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
