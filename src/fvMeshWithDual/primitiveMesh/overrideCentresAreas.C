/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Class
    Foam::primitiveMesh

Description
    Additional functions for meshWithCentres to override the cell and face
    centres

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"

// ************************************************************************* //

namespace Foam
{

void Foam::primitiveMesh::overrideFaceCentres(const pointField& newCf) const
{
    if (!faceCentresPtr_) calcFaceCentresAndAreas();
    *faceCentresPtr_ = newCf;
}

void Foam::primitiveMesh::overrideCellCentres(const pointField& newCellCtrs) const
{
    if (!cellCentresPtr_) calcCellCentresAndVols();
    *cellCentresPtr_ = newCellCtrs;
}


void Foam::primitiveMesh::overrideFaceAreas(const vectorField& newAreas) const
{
    if (!faceAreasPtr_) calcFaceCentresAndAreas();
    *faceAreasPtr_ = newAreas;
}


void Foam::primitiveMesh::overrideCellVols(const scalarList& newV) const
{
    if (!cellVolumesPtr_) calcCellCentresAndVols();
    *cellVolumesPtr_ = newV;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
