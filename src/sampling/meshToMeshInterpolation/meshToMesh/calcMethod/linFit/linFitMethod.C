/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "linFitMethod.H"
#include "pointIndexHit.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"
#include "polyFit.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linFitMethod, 0);
    addToRunTimeSelectionTable
    (
        meshToMeshMethod,
        linFitMethod,
        components
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::linFitMethod::calculateAddressing
(
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    List<DynamicList<label> > srcToTgtAddr(src_.nCells());
    List<DynamicList<scalar> > srcToTgtWght(src_.nCells());

    List<DynamicList<label> > tgtToSrcAddr(tgt_.nCells());
    List<DynamicList<scalar> > tgtToSrcWght(tgt_.nCells());

    do
    {
        // find nearest tgt cell
        findNearestCell(src_, tgt_, srcCellI, tgtCellI);

        // store src/tgt cell pair
        srcToTgtAddr[srcCellI].append(tgtCellI);
        tgtToSrcAddr[tgtCellI].append(srcCellI);

        // mark source cell srcCellI and tgtCellI as matched
        mapFlag[srcCellI] = false;

        // find new source cell
        setNextNearestCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            mapFlag,
            srcCellIDs
        );
    }
    while (srcCellI >= 0);
    
    // for the case of multiple source cells per target cell, select the
    // nearest source cell only and discard the others
    const vectorField& srcCc = src_.cellCentres();
    const vectorField& tgtCc = tgt_.cellCentres();

    forAll(tgtToSrcAddr, targetCellI)
    {
        if (tgtToSrcAddr[targetCellI].size() > 1)
        {
            const vector& tgtC = tgtCc[tgtCellI];

            DynamicList<label>& srcCells = tgtToSrcAddr[targetCellI];

            label srcCellI = srcCells[0];
            scalar d = magSqr(tgtC - srcCc[srcCellI]);

            for (label i = 1; i < srcCells.size(); i++)
            {
                label srcI = srcCells[i];
                scalar dNew = magSqr(tgtC - srcCc[srcI]);
                if (dNew < d)
                {
                    d = dNew;
                    srcCellI = srcI;
                }
            }

            srcCells.clear();
            srcCells.append(srcCellI);
        }
    }

    // If there are more target cells than source cells, some target cells
    // might not yet be mapped
    forAll(tgtToSrcAddr, tgtCellI)
    {
        if (tgtToSrcAddr[tgtCellI].empty())
        {
            label srcCellI = findMappedSrcCell(tgtCellI, tgtToSrcAddr);

            findNearestCell(tgt_, src_, tgtCellI, srcCellI);

            tgtToSrcAddr[tgtCellI].append(srcCellI);
        }
    }
    
    // Add neighbour cells to addressing
    forAll(srcToTgtAddr, i)
    {
        srcToTgtAddr[i].append(tgt().cellCells()[srcToTgtAddr[i][0]]);
    }
    Info << "tgtToSrcAddr size = " << flush << tgtToSrcAddr.size() << endl;
    forAll(tgtToSrcAddr, i)
    {
        tgtToSrcAddr[i].append(src().cellCells()[tgtToSrcAddr[i][0]]);
    }
    
    // Calculate the source to target weights based on a linear fit
    forAll(srcToTgtWght, i)
    {
        srcToTgtWght[i] = scalarList(srcToTgtAddr.size());
        
        polyFit<oONE> linFit;
        linFit.calcWeights
        (
            srcToTgtWght[i], tgt(), srcToTgtAddr[i], src().cellCentres()[i]
        );
    }

    // Calculate the target to source weights based on a quadratic fit
    forAll(tgtToSrcWght, i)
    {
        tgtToSrcWght[i] = scalarList(tgtToSrcAddr.size());
        
        polyFit<oONE> linFit;
        linFit.calcWeights
        (
            tgtToSrcWght[i], src(), tgtToSrcAddr[i], tgt().cellCentres()[i]
        );
    }

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellAddr[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellAddr[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght[i].transfer(tgtToSrcWght[i]);
    }

//    const label celli = 909;
//    Info << "srcToTgt Addressing for cell " << celli  << " "
//         << srcToTgtCellAddr[celli]
//         << "\nsrcToTgt Weights for cell " << celli  << " "
//         << srcToTgtCellWght[celli] << endl;

//    Info << "tgtToSrc Addressing for cell " << celli  << " "
//         << tgtToSrcCellAddr[celli]
//         << "\ntgtToSrc Weights for cell " << celli  << " "
//         << tgtToSrcCellWght[celli] << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linFitMethod::linFitMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    mapNearestMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::linFitMethod::~linFitMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linFitMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        calculateAddressing
        (
            srcToTgtAddr,
            srcToTgtWght,
            tgtToSrcAddr,
            tgtToSrcWght,
            srcSeedI,
            tgtSeedI,
            srcCellIDs,
            mapFlag,
            startSeedI
        );
    }
    else
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }
}


// ************************************************************************* //
