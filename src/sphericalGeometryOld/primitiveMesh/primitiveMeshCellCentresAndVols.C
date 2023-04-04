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

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.
    Cell centres and volumes for spherical geometry

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
const scalar RADIUS_TOL = 1e-5;
const scalar ANGLE_TOL = 1e-5;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    //if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void Foam::primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = Zero;
    cellVols = 0.0;

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Calculate exact cell volume and centre

    scalarField r1(nCells(), -1);
    scalarField r2(nCells(), -1);
    
    forAll (faces(), faceI)
    {
        if (!isRadialFace(fCtrs[faceI], fAreas[faceI]))
        {
            cellVols[own[faceI]] += fCtrs[faceI] & fAreas[faceI];
            cellCtrs[own[faceI]] += fCtrs[faceI];
            if (r1[own[faceI]] < 0)
            {
                r1[own[faceI]] = mag(fCtrs[faceI]);
            }
            else
            {
                r2[own[faceI]] = mag(fCtrs[faceI]);
            }
            
            if (faceI < nInternalFaces())
            {
                cellVols[nei[faceI]] += fCtrs[faceI] & fAreas[faceI];
                cellCtrs[nei[faceI]] += fCtrs[faceI];
                r2[nei[faceI]] = mag(fCtrs[faceI]);
            }
        }
    }
    cellVols *= 1./3.;
    cellCtrs *= 0.5*(r1+r2)/mag(cellCtrs);
}


bool Foam::primitiveMesh::isRadialFace
(
    const vector Cf,
    const vector Sf
) const
{
    return (sqr(Sf & Cf)/(magSqr(Sf)*magSqr(Cf))) < SMALL;
}


void Foam::primitiveMesh::overrideCellCentres(const pointField& newCellCtrs)
{
    Info << "Overriding spherical cell centres\n" << endl;

    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }
    vectorField& cellCtrs = *cellCentresPtr_;
    
    forAll(cellCtrs, cellI)
    {
        cellCtrs[cellI] = unitVector(newCellCtrs[cellI])*mag(cellCtrs[cellI]);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}



// ************************************************************************* //
