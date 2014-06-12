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

#include "CFCFaceToFaceStencil.H"
#include "syncTools.H"
#include "emptyPolyPatch.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per face the neighbour data (= face or boundary in global
// numbering). First element is always face itself!
void Foam::CFCFaceToFaceStencil::calcCellStencil
(
    labelListList& elements,
    List<labelPairList>& transformedElements
) const
{
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();


    // Non-empty boundary faces
    boolList validFace(mesh().nFaces(), true);

    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<emptyPolyPatch>(pp))
        {
            label faceI = pp.start();
            forAll(pp, i)
            {
                validFace[faceI++] = false;
            }
        }
    }


    // Calculate faces of coupled neighbour (in global numbering)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<labelPairList> neiGlobal(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        //if (!isA<emptyPolyPatch>(pp))
        if (pp.coupled())
        {
            const labelUList& faceCells = pp.faceCells();
            forAll(faceCells, i)
            {
                label faceI = pp.start()+i;
                label bFaceI = faceI-mesh().nInternalFaces();

                const cell& cFaces = mesh().cells()[faceCells[i]];

                label n = 0;
                forAll(cFaces, cFaceI)
                {
                    label fI = cFaces[cFaceI];
                    if (fI != faceI && validFace[fI])
                    {
                        n++;
                    }
                }
                neiGlobal[bFaceI].setSize(n);
                n = 0;
                forAll(cFaces, cFaceI)
                {
                    label fI = cFaces[cFaceI];
                    if (fI != faceI && validFace[fI])
                    {
                        neiGlobal[bFaceI][n++] =
                            globalIndexAndTransform::encode
                            (
                                Pstream::myProcNo(),
                                fI,
                                globalTransforms().addToTransformIndex
                                (
                                    globalTransforms().nullTransformIndex(),
                                    patchI,
                                    true           // patchI is sending side
                                )
                            );
                    }
                }
            }
        }
    }

    syncTools::syncBoundaryFaceList
    (
        mesh(),
        neiGlobal,
        eqOp<labelPairList>(),
        dummyTransform()
    );


    // Determine faces of owner and neighbour in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    elements.setSize(mesh().nFaces());

    for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
    {
        label own = mesh().faceOwner()[faceI];
        const cell& ownFaces = mesh().cells()[own];

        label nei = mesh().faceNeighbour()[faceI];
        const cell& neiFaces = mesh().cells()[nei];

        label n = 0;
        forAll(ownFaces, ownI)
        {
            if (ownFaces[ownI] != faceI && validFace[ownFaces[ownI]])
            {
                n++;
            }
        }
        forAll(neiFaces, cFaceI)
        {
            if (neiFaces[cFaceI] != faceI && validFace[neiFaces[cFaceI]])
            {
                n++;
            }
        }
        labelList& fFaces = elements[faceI];
        fFaces.setSize(n+1);
        n = 0;
        // include this face first
        fFaces[n++] = globalNumbering().toGlobal(faceI);
        forAll(ownFaces, ownI)
        {
            if (ownFaces[ownI] != faceI && validFace[ownFaces[ownI]])
            {
                fFaces[n++] = globalNumbering().toGlobal(ownFaces[ownI]);
            }
        }
        forAll(neiFaces, neiI)
        {
            if (neiFaces[neiI] != faceI && validFace[neiFaces[neiI]])
            {
                fFaces[n++] = globalNumbering().toGlobal(neiFaces[neiI]);
            }
        }
    }


    // Boundary faces

    transformedElements.setSize(mesh().nFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!isA<emptyPolyPatch>(pp))
        {
            const labelUList& faceCells = pp.faceCells();

            forAll(faceCells, i)
            {
                label faceI = pp.start()+i;
                
                const cell& ownFaces = mesh().cells()[faceCells[i]];

                DynamicList<label> untrafoFaces(ownFaces.size());
                DynamicList<labelPair> trafoFaces(ownFaces.size());

                // include this face first
                untrafoFaces.append(globalNumbering().toGlobal(faceI));

                // Collect owner side faces
                forAll(ownFaces, ownI)
                {
                    if (ownFaces[ownI] != faceI && validFace[ownFaces[ownI]])
                    {
                        untrafoFaces.append
                        (
                            globalNumbering().toGlobal(ownFaces[ownI])
                        );
                    }
                }


                if (pp.coupled())
                {
                    // Collect coupled neighbour side faces into untransformed
                    // and transformed bin.
                    const labelPairList& nbrInfo =
                        neiGlobal[faceI-mesh().nInternalFaces()];
                    forAll(nbrInfo, j)
                    {
                        // Extract transformation, index, processor from
                        // neighbour info.
                        const labelPair& info = nbrInfo[j];
                        label transform =
                            globalIndexAndTransform::transformIndex
                            (
                                info
                            );

                        if
                        (
                            transform
                         == globalTransforms().nullTransformIndex()
                        )
                        {
                            label procI =
                                globalIndexAndTransform::processor(info);
                            label index =
                                globalIndexAndTransform::index(info);
                            untrafoFaces.append
                            (
                                globalNumbering().toGlobal(procI, index)
                            );
                        }
                        else
                        {
                            trafoFaces.append(info);
                        }
                    }
                }

                elements[faceI].transfer(untrafoFaces);
                transformedElements[faceI].transfer(trafoFaces);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CFCFaceToFaceStencil::CFCFaceToFaceStencil(const polyMesh& mesh)
:
    faceToFaceStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    calcCellStencil(elements_, transformedElements_);
}


// ************************************************************************* //
