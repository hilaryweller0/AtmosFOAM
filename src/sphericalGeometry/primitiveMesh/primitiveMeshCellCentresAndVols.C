/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    primitiveMesh

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.
    Cell centres and volumes for spherical geometry

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const scalar RADIUS_TOL = 1e-5;
const scalar ANGLE_TOL = 1e-5;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcCellCentresAndVols() const
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
        FatalErrorIn("primitiveMesh::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << exit(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);
    //Info << "cellVols = " << cellVols << endl;

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = vector::zero;
    cellVols = 0.0;
    
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // next calculate exact cell volume and centre

    scalarField r1(nCells(), -1);
    scalarField r2(nCells(), -1);
    
    forAll (faces(), faceI)
    {
        scalar cos2CfAf = sqr(fAreas[faceI] & fCtrs[faceI])
                         /(magSqr(fAreas[faceI])*magSqr(fCtrs[faceI]));
        
        // determine if the face is on a sphere and hence contributes
        if(mag(cos2CfAf) >= SMALL)
        {
            cellVols[own[faceI]] += fCtrs[faceI] & fAreas[faceI];
            cellCtrs[own[faceI]] += fCtrs[faceI];
            if (r1[own[faceI]] < 0) r1 = mag(fCtrs[faceI]);
            else r2[own[faceI]] = mag(fCtrs[faceI]);
            
            if (faceI < nInternalFaces())
            {
                cellVols[nei[faceI]] += fCtrs[faceI] & fAreas[faceI];
                cellCtrs[nei[faceI]] += fCtrs[faceI];
                r2[nei[faceI]] = mag(fCtrs[faceI]);
                
                //Info << "Face " << faceI << " cos2CfAf = " << cos2CfAf << endl;
            }
        }
    }
    cellVols *= 1./3.;
    cellCtrs *= sqrt((sqr(r1) + r1*r2 + sqr(r2))/3.)/mag(cellCtrs);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
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

const scalarField& primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
