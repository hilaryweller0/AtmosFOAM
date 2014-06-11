/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "meshWithDual.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshWithDual::calcSphericalVolGeom() const
{
    // calculate faceCentres, cellCentres. faceAreas, cellVolumes, deltaCoeffs
    // based on readCentres_
    
    // First calculate faceCentres and faceAreas
    pointField fCtr = faceCentres();
    vectorField fArea = faceAreas();
    forAll(fCtr, faci)
    {
        const face& f = faces()[faci];
    
        // Check if this is a vertical face or a face on a sphere
        bool faceOnSphere = (faceToPatchFace()[faci] != 0);
        
        // Points of the face
        pointField p(f.size());
        scalarList radii(f.size());
        forAll(f, ip)
        {
            p[ip] = unitVector(points()[f[ip]]);
            radii[ip] = mag(points()[f[ip]]);
        }
        
        // Set the area for a vertical face
        if (!faceOnSphere)
        {
            // find the outer and inner radii of the face and the two rhat
            scalar a = -1, b = -1, c = -1, d = -1;
            Pair<vector> rhats;
            label rhati = 0;
            label iprev = p.size()-1;
            for(label ip = 0; ip < p.size(); ip++)
            {
                // Check if edge is radial
                if (mag(p[ip]-p[iprev]) < SMALL)
                {
                    rhats[rhati] = p[ip];
                    if (rhati == 0)
                    {
                        a = min(radii[ip], radii[iprev]);
                        b = max(radii[ip], radii[iprev]);
                    }
                    else
                    {
                        d = min(radii[ip], radii[iprev]);
                        c = max(radii[ip], radii[iprev]);
                    }
                    rhati++;
                }
                iprev = ip;
            }
            
            // rhat for the face centre
            const vector rhat = readCentres_ ? unitVector(fCtr[faci])
                                             : unitVector(rhats[0]+rhats[1]);
            // Set the face centre and area
            fCtr[faci] = 0.25*(a+b+c+d)*rhat;
            const vector idir = unitVector(rhat ^ (rhats[0] - rhats[1]));
            const scalar A = arcLength(rhats[0], rhats[1])/6.*
            (
                sqr(c) + b*c + sqr(b) - (sqr(a) + a*d + sqr(d))
            );
            fArea[faci] = A*idir*sign(fArea[faci] & idir);
        }
        else // face on the sphere
        {
            // find total solid angle using an initial approximate centre
            point fCentre = unitVector(fCtr[faci]);
            
            // next find exact centre and the total solid angle, omega
            scalar sumR = 0;
            scalar omega = 0;

            // Circulate around points on the face to calculate the angle
            label iprev = p.size()-1;
            for(label ip = 0; ip < p.size(); ip++)
            {
                scalar r = radii[ip] + radii[iprev];
                scalar a = sphTriSolidAngle(p[ip], p[iprev], fCentre);
                omega += a;
                sumR += r*a;
                iprev = ip;
            }
            if (mag(omega) < SMALL)
            {
                FatalErrorIn("calcSphericalVolGeom")
                    << " face " << faci << " has a solid anlge " << omega
                    << exit(FatalError);
            }
                
            sumR /= (2*omega);
            vector rhat = unitVector(fCtr[faci]);
            fCtr[faci] = rhat*sumR;
            fArea[faci] = sqr(sumR)*rhat*omega*sign(fArea[faci] & rhat);
            
        }
    }
    overrideFaceCentres(fCtr);
    overrideFaceAreas(fArea);
        
    // Cell centres and volumes
    pointField cellCtrs(nCells(), vector::zero);
    scalarList cellVols(nCells(), scalar(0));
    
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Sum of solid angles, alpha, for each cell
    scalarField alphaSum(nCells(), scalar(0));

    forAll(own, facei)
    {
        // Calculate volume of sector of the sphere
        scalar SfdotCf = fArea[facei] & fCtr[facei];
        // Accumulate onto owner and neighbour cells
        cellVols[own[facei]] += SfdotCf;
        if (facei < nInternalFaces())
        {
            cellVols[nei[facei]] -= SfdotCf;
        }
    
        // Calculate solid angle, alpha
        scalar alpha = mag(SfdotCf)/pow(magSqr(fCtr[facei]), 1.5);
    
        // Accumulate solid angles and centres for owner and neighbour cells
        cellCtrs[own[facei]] += alpha*fCtr[facei];
        alphaSum[own[facei]] += alpha;
        if (facei < nInternalFaces())
        {
            cellCtrs[nei[facei]] += alpha*fCtr[facei];
            alphaSum[nei[facei]] += alpha;
        }
    }

    const scalar third = 1./3.;
    forAll(cellVols, i) { cellVols[i] *= third; }
    cellCtrs /= alphaSum;
    
    if (readCentres_)
    {
        cellCtrs = unitVector(C().internalField())*mag(cellCtrs);
    }

    overrideCellCentres(cellCtrs);
    overrideCellVols(cellVols);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
