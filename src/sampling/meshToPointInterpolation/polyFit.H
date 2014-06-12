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

Class
    Foam::polyFit

Description
    Calculates the weights for a polynomial fit given a stencil of points
    withing a mesh and a point to interpolate onto

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef polyFit_H
#define polyFit_H

#include "approxType.H"
#include "ListOps.H"
#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyFit Declaration
\*---------------------------------------------------------------------------*/

template<unsigned Order>
class polyFit
:
    public approxType
{
private:
    // Private member function
        //- The minimum size for an dim dimensional polynomial of this order
        inline int minPolySize(const int dim) const
        {
            if (Order < 0 || Order > 2)
            {
                FatalErrorIn("Foam::polyFit::minPolySize")
                    << "Order must be 0, 1 or 2, not "
                    << Order << exit(FatalError);
            }

            return Order == 0 ? 1 :
              Order == 1 ? (dim == 1 ? 2 : dim == 2 ? 3 : dim == 3 ? 4 : 0) :
              Order == 2 ? (dim == 1 ? 3 : dim == 2 ? 6 : dim == 3 ? 10: 0) :
                       0;
        }


public:
        //- Destructor
        virtual ~polyFit(){}


        //- Member functions
        void calcWeights
        (
            scalarList& weights, const polyMesh& mesh,
            const labelList& stencil, const point& target
        ) const
        {
            const int nDims = mesh.nGeometricD();
            const int nCentral = 1;
            const scalar centralWeight = 1000;
            const int minSize = minPolySize(nDims);
            const label N = stencil.size();
            vectorField pts(N);
            forAll(pts, ip)
            {
                pts[ip] = mesh.cellCentres()[stencil[ip]] - target;
            }
            scalarList dists2 = magSqr(pts);
            weights.setSize(N);
            
            // Check for co-incidence with one point
            label closePoint = findMin(dists2);
            if (dists2[closePoint] < SMALL*max(dists2))
            {
                weights = 0;
                weights[closePoint] = 1;
            }
            else
            {
                // local co-ordinate directions
                label farPoint = findMax(dists2);
                vector idir(pts[farPoint]);
                idir /= mag(idir);
                // next find the furthest point away at 90 degrees to farPoint
                scalar jDist2 = magSqr(pts[0] ^ idir);
                label jPoint = 0;
                for(label j = 1; j < pts.size(); j++)
                {
                    if (j != farPoint)
                    {
                        scalar dist2 = magSqr(pts[j] ^ idir);
                        if (dist2 > jDist2)
                        {
                            jDist2 = dist2;
                            jPoint = j;
                        }
                    }
                }
                vector jdir = pts[jPoint] - (pts[jPoint] & idir)*idir;
                jdir /= mag(jdir);

                vector kdir = idir ^ jdir;

                // scale factor for the geometry
                scalar scale = sqrt(max(dists2));

                // Re-calculate the list of points in the local co-ordinate
                // directions scaled by the geometry scaling factor.
                forAll(pts, ip)
                {
                    const point p = pts[ip];
                    pts[ip].x() = (p & idir)/scale;
                    pts[ip].y() = nDims >= 2 ? (p & jdir)/scale : 0;
                    pts[ip].z() = nDims >= 3 ? (p & kdir)/scale : 0;
                }

                // weights for different points
                scalarList wts(N, scalar(1));
                for(label iw = 0; iw < nCentral; iw++) wts[iw] = centralWeight;

                // calculate the matrix of the polynomial components
                scalarRectangularMatrix B(N, minSize, scalar(0));

                forAll(pts, ip)
                {
                    scalar px = pts[ip].x();
                    scalar py = pts[ip].y();
                    scalar pz = pts[ip].z();

                    label is = 0;
                    B[ip][is++] = wts[ip];

                    if (Order >= 1) B[ip][is++] = wts[ip]*px;
                    if (Order >= 2) B[ip][is++] = wts[ip]*sqr(px);

                    if (nDims >= 2)
                    {
                        if (Order >= 1) B[ip][is++] = wts[ip]*py;
                        if (Order >= 2)
                        {
                            B[ip][is++] = wts[ip]*px*py;
                            B[ip][is++] = wts[ip]*sqr(py);
                        }
                    }
                    if (nDims == 3)
                    {
                        if (Order >= 1) B[ip][is++] = wts[ip]*pz;
                        if (Order >= 2)
                        {
                            B[ip][is++] = wts[ip]*px*pz;
                            B[ip][is++] = wts[ip]*py*pz;
                            B[ip][is++] = wts[ip]*sqr(pz);
                        }
                    }
                }

                SVD svd(B, SMALL);

                // Set the coefficients
                forAll(weights, i)
                {
                    weights[i] = wts[i]*svd.VSinvUt()[0][i];
                }
            }
            if (mag(sum(weights)-1) > 1e-6)
            {
                FatalErrorIn("polyFit::calcWeights")
                << " the weights should sum to 1 but weights = " << weights
                    << " and 1-sum(weights) = " << 1 - sum(weights)
                    << exit(FatalError);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
