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

//#include "volPointInterpolation.H"
//#include "volFields.H"
#include "pointConstraints.H"
#include "coupledPointPatchField.H"
//#include "cyclicPointPatch.H"

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

    const surfaceVectorField delta = mesh.delta();
    
    tmp<GeometricField<GradType, pointPatchField, pointMesh> > treconField
    (
        new GeometricField<GradType, pointPatchField, pointMesh>
        (
            IOobject("volIntegrate("+ssf.name()+')', ssf.instance(), mesh),
            pointMesh::New(mesh),
            dimensioned<GradType>("", ssf.dimensions(), GradType::zero),
            "calculated"
        )
    );
    GeometricField<GradType, pointPatchField, pointMesh>& grad = treconField.ref();
    
    // The tensor to invert to find the least squares reconstruction
    GeometricField<tensor, pointPatchField, pointMesh> invTensor
    (
        IOobject("invTensor", ssf.instance(), mesh),
        pointMesh::New(mesh),
        dimensioned<tensor>("", dimless, tensor::zero),
        "calculated"
    );
    
    // Loop over all faces and add contribution of gradeint to the points
    // and contribution to the tensor
    forAll(ssf, faceI)
    {
        // Use an extended stencil for a mesh whose dual is a triangulation
        bool extend
        (
            mesh.cells()[mesh.owner()[faceI]].size() >= 7
         && mesh.cells()[mesh.neighbour()[faceI]].size() >= 7
        );
    
        const vector& d = delta[faceI];
        const vector dhat = d/mag(d);
    
        // Find all points of the face and loop over them
        const labelList& face = mesh.faces()[faceI];
        forAll(face, ip)
        {
            label iPt = face[ip];
            grad[iPt] += ssf[faceI] * dhat;
            invTensor[iPt] += dhat*dhat;

            // Include a larger stencil for a mesh which is the dual
            // of a triangulation
            if (extend)
            {
                // Find all points of the faces of the points
                const labelList& pointFaces = mesh.pointFaces()[iPt];
                forAll(pointFaces, ipf)
                {
                    const label iff = pointFaces[ipf];
                    const labelList face2 = mesh.faces()[iff];
                    forAll(face2, ipp)
                    {
                        label iPtf = face2[ipp];
                        grad[iPtf] += ssf[faceI] * dhat;
                        invTensor[iPtf] += dhat * dhat;
                    }
                }
            }
        }
        
    }

    // And all boundary faces
    forAll(ssf.boundaryField(), patchI)
    {
        const scalarField& ssfb = ssf.boundaryField()[patchI];
        forAll(ssfb, faceI)
        {
            const vector& d = delta.boundaryField()[patchI][faceI];
            const vector dhat = d/mag(d);
            
            const labelList& face = mesh.boundaryMesh()[patchI][faceI];
            forAll(face, ip)
            {
                label iPt = face[ip];
                grad[iPt] += ssfb[faceI]*dhat;
                invTensor[iPt] += dhat*dhat;
            }
        }
    }

    // Apply coupled boundary conditions to invTensor and to grad
    applyCoupledBCs(grad, mesh);
    applyCoupledBCs(invTensor, mesh);

    // Calculate the gradient on the points
    grad = inv(invTensor) & grad;

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
         tf(fvc::faceToPointReconstruct(tssf.ref()));
    tssf.clear();
    return tf;
}


template<class Type>
void applyCoupledBCs
(
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const fvMesh& mesh
)
{
    // Apply boundary conditions to pf
    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pfi = pf.ref();

    pointConstraints::syncUntransformedData(mesh, pfi, plusEqOp<Type>());
    addSeparated(pf);
    pushUntransformedData(pf, mesh);
}


template<class Type>
void addSeparated(GeometricField<Type, pointPatchField, pointMesh>& pf)
{
    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pfi = pf.ref();

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& pfbf = pf.boundaryFieldRef();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).initSwapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }

    // Block for any outstanding requests
    Pstream::waitRequests();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).swapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }
}


template<class Type>
void pushUntransformedData
(
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const fvMesh& mesh
)
{
    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pointData = pf.ref();

    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
