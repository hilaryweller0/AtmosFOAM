/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "approxType.H"

// * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

void Foam::approxType::findDirs
(
    List<vector>& dir, const pointField& points,
    const labelList& stencil, const point& target, const int nDims
) const
{
    dir.setSize(3);

    // Points in local coordinate system
    const label N = stencil.size();
    vectorField pts(N);
    forAll(pts, ip)
    {
        pts[ip] = points[stencil[ip]] - target;
    }
    scalarList dists2 = magSqr(pts);

    // local co-ordinate directions
    label farPoint = findMax(dists2);
    dir[0] = pts[farPoint];
    dir[0] /= mag(dir[0]);
    // next find the furthest point away at 90 degrees to farPoint
    scalar jDist2 = magSqr(pts[0] ^ dir[0]);
    label jPoint = 0;
    for(label j = 1; j < pts.size(); j++)
    {
        if (j != farPoint)
        {
            scalar dist2 = magSqr(pts[j] ^ dir[0]);
            if (dist2 > jDist2)
            {
                jDist2 = dist2;
                jPoint = j;
            }
        }
    }
    dir[1] = pts[jPoint] - (pts[jPoint] & dir[0])*dir[0];
    dir[1] /= mag(dir[1]);

    dir[2] = dir[0] ^ dir[1];
}


// ************************************************************************* //
