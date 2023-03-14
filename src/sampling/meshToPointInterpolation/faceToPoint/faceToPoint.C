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

#include "faceToPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class approxType>
Foam::faceToPoint<approxType>::faceToPoint
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    weights_(mesh.nPoints())
{
    labelList pointFaces(3);
    forAll(weights_, ip)
    {
        const point& p = mesh.points()[ip];
        //const labelList& pointFaces = mesh.pointFaces()[ip];
        pointFaces[0] = mesh.pointFaces()[ip][0];
        pointFaces[1] = mesh.pointFaces()[ip][1];
        pointFaces[2] = mesh.pointFaces()[ip][2];
        weights_[ip].setSize(pointFaces.size());

        // Create the weights (depending on the approxType)
        approxType thisApproxType;
        thisApproxType.calcWeights
        (
            weights_[ip],
            mesh_.faceCentres(),
            pointFaces,
            p,
            2
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class approxType>
Foam::faceToPoint<approxType>::~faceToPoint()
{}


// ************************************************************************* //
