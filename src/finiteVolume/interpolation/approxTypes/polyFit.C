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

#include "polyFit.H"

// * * * * * * * * * * * Private Member Function * * * * * * * * * * * * * * //

template<Foam::orderType Order>
int Foam::polyFit<Order>::minPolySize
(
    const int dim,
    const orderType order
) const
{
    return order == oZERO ? 1 :
      order == oONE ?   (dim == 1 ? 2 : dim == 2 ? 3  : dim == 3 ? 4 : 0) :
      order == oTWO ?   (dim == 1 ? 3 : dim == 2 ? 6  : dim == 3 ? 10: 0) :
      order == oTHREE ?   (dim == 1 ? 4 : dim == 2 ? 10 : dim == 3 ? 20: 0) :
      order == oONEPLUS ? (dim == 1 ? 3 : dim == 2 ? 5  : dim == 3 ? 7: 0) :
      order == oTWOPLUS ? (dim == 1 ? 3 : dim == 2 ? 8  : dim == 3 ? 15: 0) :
      order == oTWOPLUSPLUS ? (dim == 1 ? 4 : dim == 2 ? 9  : dim == 3 ? 18: 0) :
      order == oTHREEPLUS ? (dim == 1 ? 5 : dim == 2 ? 13 : dim == 3 ? 26: 0) :
             0;
}


// * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

template<Foam::orderType Order>
int Foam::polyFit<Order>::calcWeights
(
    scalarList& weights, const polyMesh& mesh,
    const labelList& stencil, const point& target, const int centralCell
) const
{
    const int nDims = mesh.nSolutionD();
    List<vector> dir(3);
    findDirs(dir, mesh.cellCentres(), stencil, target, nDims);
    return calcWeights(weights, mesh, stencil, target, dir, centralCell);
}


template<Foam::orderType Order>
int Foam::polyFit<Order>::calcWeights
(
    scalarList& weights, const polyMesh& mesh,
    const labelList& stencil, const point& target,
    const List<vector>& dir, const int centralCell
) const
{
    const int nDims = mesh.nSolutionD();
    return calcWeights
    (
        weights, mesh.cellCentres(), stencil, target, dir, nDims,
        centralCell
    );
}


template<Foam::orderType Order>
int Foam::polyFit<Order>::calcWeights
(
    scalarList& weights, const pointField& points,
    const labelList& stencil, const point& target,
    const int nDims, const int centralCell
) const
{
    List<vector> dir(3);
    findDirs(dir, points, stencil, target, nDims);
    return calcWeights
    (
        weights, points, stencil, target, dir, nDims, centralCell
    );
}


template<Foam::orderType Order>
int Foam::polyFit<Order>::calcWeights
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


template<Foam::orderType Order>
int Foam::polyFit<Order>::calcWeights
(
    scalarList& weights, const pointField& points,
    const labelList& stencil, const point& target,
    const List<vector>& dir, const int nDims, const int centralCell
) const
{
    const int nCentral = 1;
    const scalar centralWeight = 1000;
    int minSize = minPolySize(nDims);
    int ord = Order;
    const label N = stencil.size();
    vectorField pts(N);
    forAll(pts, ip)
    {
        pts[ip] = points[stencil[ip]] - target;
    }
    scalarList dists2 = magSqr(pts);
    weights.setSize(N);
    
    while (N < minSize)
    {
//        WarningIn("Foam::polyFit<Order>::calcWeights") << "nDims = " << nDims
//            << " Order = " << Order << " min stencil size = " << minSize
//            << " stencilSize = " << N << " reducing order" << endl;
        minSize = minPolySize(nDims, static_cast<orderType>(--ord));
    }

    // scale factor for the geometry
    scalar scale = sqrt(max(dists2));

    // weights for goodness of fit for different points
    scalarList wts(N, scalar(1));
    for(label iw = 0; iw < nCentral; iw++)
    {
        wts[iw] = centralWeight;
    }
    //if (centralCell != -1) wts[centralCell] = centralWeight;

    // Re-calculate the list of points in the local co-ordinate
    // directions scaled by the geometry scaling factor.
    forAll(pts, ip)
    {
        const point p = pts[ip];
        pts[ip].x() = (p & dir[0])/scale;
        pts[ip].y() = nDims >= 2 ? (p & dir[1])/scale : 0;
        pts[ip].z() = nDims >= 3 ? (p & dir[2])/scale : 0;
    }

    // Check for co-incidence with one point
    label closePoint = findMin(dists2);
    if (dists2[closePoint] < SMALL*max(dists2))
    {
        weights = 0;
        weights[closePoint] = 1;
    }
    else for
    (
        bool goodFit = false;
        !goodFit;
        minSize = minPolySize(nDims, static_cast<orderType>(--ord))
    )
    {
        orderType order = static_cast<orderType>(ord);
    
        // calculate the matrix of the polynomial components
        scalarRectangularMatrix B(N, minSize, scalar(0));

        forAll(pts, ip)
        {
            scalar px = pts[ip].x();
            scalar py = pts[ip].y();
            scalar pz = pts[ip].z();

            label is = 0;
            B[ip][is++] = wts[ip];

            if (order >= oONE)   B[ip][is++] = wts[ip]*px;
            if (order >= oONEPLUS) B[ip][is++] = wts[ip]*sqr(px);
            if (order >= oTWOPLUSPLUS)   B[ip][is++] = wts[ip]*pow(px,3);
            if (order >= oTHREEPLUS) B[ip][is++] = wts[ip]*pow(px,4);

            if (nDims >= 2)
            {
                if (order >= oONE) B[ip][is++] = wts[ip]*py;
                if (order >= oONEPLUS) B[ip][is++] = wts[ip]*sqr(py);
                if (order >= oTWO)
                {
                    B[ip][is++] = wts[ip]*px*py;
                }
                if (order >= oTWOPLUS)
                {
                    B[ip][is++] = wts[ip]*sqr(px)*py;
                    B[ip][is++] = wts[ip]*px*sqr(py);
                }
                if (order >= oTHREE)
                {
                    B[ip][is++] = wts[ip]*pow(py,3);
                }
                if (order >= oTHREEPLUS)
                {
                    B[ip][is++] = wts[ip]*pow(py,4);
                    B[ip][is++] = wts[ip]*sqr(px)*sqr(py);
                }
            }
            if (nDims == 3)
            {
                if (order >= oONE) B[ip][is++] = wts[ip]*pz;
                if (order >= oTWO)
                {
                    B[ip][is++] = wts[ip]*px*pz;
                    B[ip][is++] = wts[ip]*py*pz;
                    B[ip][is++] = wts[ip]*sqr(pz);
                }
                if (order >= oTWOPLUSPLUS)
                {
                    B[ip][is++] = wts[ip]*px*py*pz;
                    B[ip][is++] = wts[ip]*px*sqr(py);
                    B[ip][is++] = wts[ip]*px*sqr(pz);
                    B[ip][is++] = wts[ip]*sqr(px)*py;
                    B[ip][is++] = wts[ip]*py*sqr(pz);
                    B[ip][is++] = wts[ip]*px*sqr(pz);
                    B[ip][is++] = wts[ip]*py*sqr(pz);
                }
                if (order >= oTHREE)
                {
                    B[ip][is++] = wts[ip]*pow(pz,3);
                }
                if (order >= oTHREEPLUS)
                {
                    B[ip][is++] = wts[ip]*pow(pz,4);
                    B[ip][is++] = wts[ip]*sqr(px)*sqr(pz);
                    B[ip][is++] = wts[ip]*sqr(py)*sqr(pz);
                }
            }
        }
        SVD svd(B, SMALL);

        // Set the coefficients
        forAll(weights, i)
        {
            weights[i] = wts[i]*svd.VSinvUt()[0][i];
        }

        // Reasons for failure of the fit and hence to drop to a lower order
        if
        (
            mag(sum(weights)-1) > 1e-6
         || max(weights) > weights[0]
         || weights[0] > 1.5
         || weights[0] <= 0 || weights[centralCell] <= 0
        )
        {
//            Info << "ord = " << ord << " singular values = " << svd.S() << endl;
//            Info << "nDims = " << nDims << " ord = " << ord
//                 << " min stencil size = " << minSize << " stencilSize = "
//                 << N << " singular values = " << svd.S()
//                << "\nB = " << B
//                 << "\nsvd.VSinvUt() = " << svd.VSinvUt() << endl;
//        
//            WarningIn("Foam::polyFit<Order>::calcWeights")
//            << " the weights should sum to 1 but weights = " << weights
//                << " and 1-sum(weights) = " << 1 - sum(weights)
//                << " reducing order"
//                << endl;
        }
        else
        {
            goodFit = true;
//            Info << "weights" << weights << "\nB = " << B
//                 << "\nsvd.VSinvUt() = " << svd.VSinvUt()<< endl;
        }
    }
    
    return ord+1;
}

// ************************************************************************* //
