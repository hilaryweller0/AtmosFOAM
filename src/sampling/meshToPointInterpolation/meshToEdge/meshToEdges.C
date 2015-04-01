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

#include "meshToEdges.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToEdges<approxType>::meshToEdges
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    weights_(mesh.nInternalEdges())
{
    forAll(weights_, ie)
    {
        const edge& e = mesh.edges()[ie];
        const point& ps = mesh.points()[e[0]];
        const point& pe = mesh.points()[e[1]];
        const point pce = 0.5*(ps + pe);
        
        const labelList& edgeCells = mesh.edgeCells()[ie];
        weights_[ie].setSize(edgeCells.size());

        // local coordinate system
        List<vector> dirs(2, vector::zero);
        vector kdir = (ps - pe)/mag(ps - pe);
        dirs[0] = mesh.C()[edgeCells[0]] - pce;
        dirs[0] -= (dirs[0] & kdir) * kdir;
        dirs[1] = kdir ^ dirs[0];

        // Create the weights (depending on the approxType)
        approxType thisApproxType;
        thisApproxType.calcWeights(weights_[ie], mesh_, edgeCells, pce, dirs);
        
//        Info << "Edge " << ie << "(" << e[0] << "," << e[1] << ") centre at "
//             << pce << " edgeCells = " << edgeCells << " weights = " << weights_[ie]
//             << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToEdges<approxType>::~meshToEdges()
{}


// ************************************************************************* //
