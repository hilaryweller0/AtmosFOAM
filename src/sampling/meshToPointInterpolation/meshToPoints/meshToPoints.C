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

#include "meshToPoints.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPoints<approxType>::meshToPoints
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    weights_(mesh.nPoints())
{
    forAll(weights_, ip)
    {
        const point& p = mesh.points()[ip];
        const labelList& pointCells = mesh.pointCells()[ip];
        weights_[ip].setSize(pointCells.size());

        // Create the weights (depending on the approxType)
        approxType thisApproxType;
        thisApproxType.calcWeights(weights_[ip], mesh_, pointCells, p);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPoints<approxType>::~meshToPoints()
{}


// ************************************************************************* //
