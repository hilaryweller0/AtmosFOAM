/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2015-11-17 AtmosFOAM, Hilary Weller, University of Reading added support for
 cyclic boundaries
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

#include "cellToFaceExtStencil.H"
#include "emptyPolyPatch.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFaceExtStencil::merge
(
    const labelPairList& lst1,
    labelPairList& lst2
)
{
    label nAdd = 0;
    forAll(lst1, i)
    {
        if (findIndex(lst2, lst1[i]) == -1)
        {
            nAdd++;
        }
    }

    label sz = lst2.size();
    lst2.setSize(sz+nAdd);

    forAll(lst1, i)
    {
        if (findIndex(lst2, lst1[i]) == -1)
        {
            lst2[sz++] = lst1[i];
        }
    }
}


void Foam::cellToFaceExtStencil::merge
(
    const label faceI,
    const labelList& ownCCells,
    const labelList& neiCCells,

    labelList& faceStencil,
    labelHashSet& workSet
)
{
    workSet.clear();

    label globalOwn = ownCCells[0];
    // Insert cellCells
    forAll(ownCCells, i)
    {
        workSet.insert(ownCCells[i]);
    }

    label globalNei = -1;
    if (neiCCells.size())
    {
        globalNei = neiCCells[0];
        // Insert cellCells
        forAll(neiCCells, i)
        {
            workSet.insert(neiCCells[i]);
        }
    }

    // Guarantee owner first, neighbour second.
    faceStencil.setSize(workSet.size());
    label n = 0;
    faceStencil[n++] = globalOwn;
    if (globalNei != -1)
    {
        faceStencil[n++] = globalNei;
    }
    forAllConstIter(labelHashSet, workSet, iter)
    {
        if (iter.key() != globalOwn && iter.key() != globalNei)
        {
            faceStencil[n++] = iter.key();
        }
    }
}


// Calculates per face a list of global cell/face indices.
void Foam::cellToFaceExtStencil::calcFaceStencil
(
    const labelListList& untrafoCellCells,
    const List<labelPairList>& trafoCellCells,
    labelListList& faceStencil,
    List<labelPairList>& trafoFaceStencil
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    // Combine (cell)stencils on either side of face. Stencil can
    // either be without transformation (untrafoCellCells) or with
    // (trafoCellCells). In case of the face being on a coupled patch
    // this can additionally add a transformation.


    // Determine neighbouring cell stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: not adding transformation here; doing that below. Could probably
    // write combineOp to add a transformation.

    labelListList neiUntrafoCellCells(nBnd);
    List<labelPairList> neiTrafoCellCells(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                label bFaceI = faceI-mesh_.nInternalFaces();
                neiUntrafoCellCells[bFaceI] = untrafoCellCells[own[faceI]];
                neiTrafoCellCells[bFaceI] = trafoCellCells[own[faceI]];
                faceI++;
            }
        }
    }
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiUntrafoCellCells,
        eqOp<labelList>(),
        dummyTransform()
    );
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiTrafoCellCells,
        eqOp<labelPairList>(),
        dummyTransform()
    );



    // Construct stencil in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    faceStencil.setSize(mesh_.nFaces());
    trafoFaceStencil.setSize(mesh_.nFaces());

    labelHashSet faceStencilSet;

    // On internal faces just merge stencils from both sides

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        // Merge owner and neighbour stencils, filter out duplicates
        // and insert in order (owner, neighbour, rest)
        merge
        (
            faceI,
            untrafoCellCells[own[faceI]],
            untrafoCellCells[nei[faceI]],

            faceStencil[faceI],
            faceStencilSet
        );

        // Do stencils with transformations
        trafoFaceStencil[faceI] = trafoCellCells[own[faceI]];
        merge(trafoCellCells[nei[faceI]], trafoFaceStencil[faceI]);
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            // Potentially add a transformation
            const labelPair& transSign =
                globalTransforms().patchTransformSign()[patchI];

            if (transSign.first() == globalTransforms().nullTransformIndex())
            {
                // The easy case: no transform
                forAll(pp, i)
                {
                    label bfaceI = faceI-mesh_.nInternalFaces();

                    // Do stencil without transformations
                    merge
                    (
                        faceI,
                        untrafoCellCells[own[faceI]],
                        neiUntrafoCellCells[bfaceI],

                        faceStencil[faceI],
                        faceStencilSet
                    );

                    // Do stencil with transformations
                    trafoFaceStencil[faceI] = trafoCellCells[own[faceI]];
                    merge(neiTrafoCellCells[bfaceI], trafoFaceStencil[faceI]);

                    faceI++;
                }
            }
            else
            {
                // Add the transformation
                forAll(pp, i)
                {
                    label bFaceI = faceI-mesh_.nInternalFaces();

                    faceStencil[faceI] = untrafoCellCells[own[faceI]];
                    trafoFaceStencil[faceI] = trafoCellCells[own[faceI]];

                    // Neighbour stencils gets transform added

                    // 1. Untransformed part of neighbour stencil
                    {
                        const labelList& neiStencil =
                            neiUntrafoCellCells[bFaceI];

                        // Add the transformation to the neiStencil
                        labelPairList transformed(neiStencil.size());
                        forAll(neiStencil, j)
                        {
                            const label elem = neiStencil[i];
                            // Get components of neighbouring data
                            label neiProc =
                                globalNumbering().whichProcID(elem);
                            label neiIndex =
                                globalNumbering().toLocal(neiProc, elem);

                            transformed[j] = globalIndexAndTransform::encode
                            (
                                neiProc,
                                neiIndex,
                                globalTransforms().addToTransformIndex
                                (
                                    globalTransforms().nullTransformIndex(),
                                    patchI,
                                    false   // patchI is the receiving patch
                                )
                            );
                        }
                        merge(transformed, trafoFaceStencil[faceI]);
                    }

                    // 2. Transformed part of neighbour stencil
                    {
                        const labelPairList& neiStencil =
                            neiTrafoCellCells[bFaceI];

                        // Add the transformation to the neiStencil
                        labelPairList transformed(neiStencil.size());
                        forAll(neiStencil, j)
                        {
                            const labelPair& elem = neiStencil[i];
                            // Get components of neighbouring data
                            label neiProc =
                                globalIndexAndTransform::processor(elem);
                            label neiIndex =
                                globalIndexAndTransform::index(elem);
                            label neiTransform =
                                globalIndexAndTransform::transformIndex(elem);

                            transformed[j] = globalIndexAndTransform::encode
                            (
                                neiProc,
                                neiIndex,
                                globalTransforms().addToTransformIndex
                                (
                                    neiTransform,
                                    patchI,
                                    false   // patchI is the receiving patch
                                )
                            );
                        }
                        merge(transformed, trafoFaceStencil[faceI]);
                    }

                    faceI++;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            // Just pick up owner side stencils.

            forAll(pp, i)
            {
                faceStencil[faceI] = untrafoCellCells[own[faceI]];
                trafoFaceStencil[faceI] = trafoCellCells[own[faceI]];

                faceI++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFaceExtStencil::cellToFaceExtStencil(const polyMesh& mesh)
:
    mesh_(mesh),
    globalNumbering_(mesh_.nCells()+mesh_.nFaces()-mesh_.nInternalFaces()),
    globalTransforms_(mesh_)
{}


// ************************************************************************* //
