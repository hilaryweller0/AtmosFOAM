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

#include "rbfFit.H"

// * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //


int Foam::rbfFit::calcWeights
(
    scalarList& weights, const polyMesh& mesh,
    const labelList& stencil, const point& target, const int centralCell
) const
{
    const int nDims = mesh.nSolutionD();
    return calcWeights
    (
        weights, mesh.cellCentres(), stencil, target, nDims, centralCell
    );
}


int Foam::rbfFit::calcWeights
(
    scalarList& weights, const pointField& points,
    const point& target, const int nDims, const int centralCell
) const
{
    labelList stencil(points.size());
    forAll(stencil, ip){ stencil[ip] = ip; }
    return calcWeights
    (
        weights, points, stencil, target, nDims, centralCell
    );
}


int Foam::rbfFit::calcWeights
(
    scalarList& weights, const pointField& points,
    const labelList& stencil, const point& target,
    const int nDims, const int centralCell
) const
{
    const label N = stencil.size();
    vectorField pts(N);
    forAll(pts, ip)
    {
        pts[ip] = points[stencil[ip]] - target;
    }
    scalarList dists2 = magSqr(pts);
    weights.setSize(N);
    
    // the RBF scale factor is the largest distance to the centre of the
    // stencil
    const scalar scale = 0.1*sqrt(max(dists2));
    
    pts /= scale;
    
    // Check for co-incidence with one point
    label closePoint = findMin(dists2);
    if (dists2[closePoint] < SMALL*max(dists2))
    {
        weights = 0;
        weights[closePoint] = 1;
    }
    else
    {
        // Matrix to invert to calculate the RBF coefficients
        scalarSquareMatrix B(N+1, N+1, scalar(0));
        // Array of rbf values
        scalarField rbfVals(N+1, scalar(0));
        rbfVals[N] = 1;
        B[N][N] = 0;

        forAll(pts, ip)
        {
            B[N][ip] = B[ip][N] = 1;
            rbfVals[ip] = rbf(mag(pts[ip]));
        
            for(label j = ip; j < N; j++)
            {
                B[ip][j] = B[j][ip] = rbf(mag(pts[ip] - pts[j]));
            }
        }

        solve(B, rbfVals);

        // Set the coefficients
        forAll(weights, i)
        {
            weights[i] = rbfVals[i];
        }
        
        if (mag(sum(weights)-1) > 1e-6)
        {
            FatalErrorIn("rbfFit::calcWeights")
            << " the weights should sum to 1 but weights = " << weights
                << " and 1-sum(weights) = " << 1 - sum(weights)
                << "\nB = " << B
                << exit(FatalError);
        }
    }
    return N;
}

// ************************************************************************* //
