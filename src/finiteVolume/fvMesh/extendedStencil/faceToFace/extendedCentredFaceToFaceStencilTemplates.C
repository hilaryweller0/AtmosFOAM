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

#include "extendedCentredFaceToFaceStencil.H"
#include "mapDistribute.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class TransformOp>
void Foam::extendedCentredFaceToFaceStencil::collectData
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    List<List<Type> >& stencilFld,
    const TransformOp& top
) const
{
    // 1. Construct flat field for internal and patch data
    List<Type> flatFld(map().constructSize(), pTraits<Type>::zero);

    // Insert my internal values
    forAll(fld, cellI)
    {
        flatFld[cellI] = fld[cellI];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& pfld = fld.boundaryField()[patchI];

        label nCompact = pfld.patch().start();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    map().distribute(fld.mesh().globalData().globalTransforms(), flatFld, top);


    // 2. Pull to stencil
    stencilFld.setSize(elements_.size());

    forAll(elements_, elemI)
    {
        const labelList& cCells = elements_[elemI];
        const labelList& trafoCCells = transformedElements_[elemI];

        stencilFld[elemI].setSize(cCells.size()+trafoCCells.size());

        label nCompact = 0;
        forAll(cCells, i)
        {
            stencilFld[elemI][nCompact++] = flatFld[cCells[i]];
        }
        forAll(trafoCCells, i)
        {
            stencilFld[elemI][nCompact++] = flatFld[trafoCCells[i]];
        }
    }
}


template<class Type>
void Foam::extendedCentredFaceToFaceStencil::collectData
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    List<List<Type> >& stencilFld
) const
{
    collectData(fld, stencilFld, mapDistribute::transform());
}


template<class Type>
void Foam::extendedCentredFaceToFaceStencil::collectPositions
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    List<List<Type> >& stencilFld
) const
{
    collectData(fld, stencilFld, mapDistribute::transformPosition());
}


template<class Type>
Foam::tmp<Foam::GeometricField<Foam::Vector<Type>, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::extendedCentredFaceToFaceStencil::weightedSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    const List<List<vector> >& stencilVectorWeights
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = fld.mesh();

    // Collect values
    List<List<Type> > stencilFld;
    collectData(fld, stencilFld);

    tmp<GeometricField<GradType, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<GradType, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<GradType>
            (
                fld.name(),
                fld.dimensions()/dimArea,
                pTraits<GradType>::zero
            )
        )
    );
    GeometricField<GradType, fvsPatchField, surfaceMesh>& sf = tsfCorr.ref();

    // Internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        const List<Type>& stField = stencilFld[faceI];
        const List<vector>& stWeight = stencilVectorWeights[faceI];

        forAll(stField, i)
        {
            sf[faceI] += stField[i]*stWeight[i];
        }
    }

    // Boundaries
    typename GeometricField<GradType, fvsPatchField, surfaceMesh>::
        Boundary& bSfCorr = sf.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<GradType>& pSfCorr = bSfCorr[patchi];

        label faceI = pSfCorr.patch().patch().start();

        forAll(pSfCorr, i)
        {
            const List<Type>& stField = stencilFld[faceI];
            const List<vector>& stWeight = stencilVectorWeights[faceI];

            forAll(stField, j)
            {
                pSfCorr[i] += stField[j]*stWeight[j];
            }

            faceI++;
        }
    }

    return tsfCorr;
}


// ************************************************************************* //
