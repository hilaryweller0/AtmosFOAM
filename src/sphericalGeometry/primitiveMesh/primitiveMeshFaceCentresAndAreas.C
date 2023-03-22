/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Description
    Calculate the face centres and areas for spherical geometry.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

    For spherical geometry - only copes with either purely radial or purely
    spherical faces - not a mixture. A spherical face is a face on the sphere
    whereas a radial face points around the sphere

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
//#include "IFstream.H"
//#include "IOmanip.H"
#include "VectorSpaceFunctions.H"

namespace Foam
{
const scalar RS_TOL = 1e-7;
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcFaceCentresAndAreas() const
{
    //if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Calculating face centres and face areas"
            << " for shallowAtmosphere spherical geometry"
            << endl;
    }

    // It is an error to attempt to recalculate faceCentres
    // if the pointer is already set
    if (faceCentresPtr_ || faceAreasPtr_ || magFaceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    magFaceAreasPtr_ = new scalarField(nFaces());
    scalarField& magfAreas = *magFaceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas, magfAreas);

    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Finished calculating face centres and face areas"
            << endl;
    }
}


void Foam::primitiveMesh::makeFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
) const
{
    // Loop through every face
    forAll (faces(), facei)
    {
        const labelList& f = faces()[facei];
        label nPoints = f.size();

        // reference radii for the face
        const scalar r1 = mag(p[f[0]]);
        const scalar r1S = sqr(r1);

        scalar r2S = r1S;
        for(label ip = 1; ip < nPoints && magSqr(r1S - r2S) < RS_TOL*r1S; ip++)
        {
            r2S = magSqr(p[f[ip]]);
        }
        const scalar r2 = sqrt(r2S);

        // Determine whether it is a radial or spherical face
        bool faceOnSphere = true;

        // for a radial face - 4 vertices with exactly 2 different lat-lon locations
        if(nPoints == 4)
        {
            const vector rhat0 = p[f[0]]/mag(p[f[0]]);
            vector rhat1 = p[f[1]]/mag(p[f[1]]);
            vector rhat2 = p[f[2]]/mag(p[f[2]]);
            vector rhat3 = p[f[3]]/mag(p[f[3]]);
            
            if (magSqr(rhat0 - rhat1) <= RS_TOL*magSqr(rhat0 - rhat2))
            {
                rhat1 = rhat2;
                rhat2 = rhat0;
            }

            // Check if there are only 2 different rhat
            const scalar distS = magSqr(rhat0 - rhat1);
            faceOnSphere = !
            (
                (
                    magSqr(rhat0 - rhat2) <= distS*RS_TOL
                 && magSqr(rhat1 - rhat3) <= distS*RS_TOL
                )
             || (
                    magSqr(rhat0 - rhat3) <= distS*RS_TOL
                 && magSqr(rhat1 - rhat2) <= distS*RS_TOL
                )
            );
            
            // Check that a radial face as 2 distinct radii
            if (!faceOnSphere && magSqr(r1S - r2S) < RS_TOL*r1S)
            {
                Info << "Face " << facei << " with points " << f[0] << " "
                     << f[1] << " " << f[2] << " " << f[3]
                     << "\nat\n"
                     << p[f[0]] << nl
                     << p[f[1]] << nl
                     << p[f[2]] << nl
                     << p[f[3]] << endl;
            
                FatalErrorIn("primitiveMesh::makeFaceCentresAndAreas")
                    << "Face " << facei
                    << " is diagnosed to be a radial face but only has one distinct radius. r1 = " << r1 << " r2 = " << r2 << " r2-r1 = " << r2-r1 << endl
                    << "rhat0 = " << rhat0 << " r = " << mag(p[f[0]]) << nl
                    << "rhat1 = " << rhat1 << " r = " << mag(p[f[1]]) << nl
                    << "rhat2 = " << rhat2 << " r = " << mag(p[f[2]]) << nl
                    << "rhat3 = " << rhat3 << " r = " << mag(p[f[3]]) << nl
                    << exit(FatalError);
            }
            
            // For a face not on the sphere (between rhat0 and rhat1) determine
            // the area and centre
            if (!faceOnSphere)
            {
                const vector rhat = (0.5*(rhat0 + rhat1))/mag(0.5*(rhat0 + rhat1));
                fCtrs[facei] = 0.5*(r1 + r2)*rhat;
                const vector idir = rhat ^ (rhat1 - rhat0);
                magfAreas[facei] = 0.5*arcLength(rhat0, rhat1)*mag(r2S - r1S);
                fAreas[facei] = magfAreas[facei]*idir/mag(idir);
                const vector dir = (p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]);
                if ((dir & idir) < 0) fAreas[facei] = -fAreas[facei];
            }
        }

        // else a spherical face
        if (faceOnSphere)
        {
            // Use an initial approximate face centre
            point fCentre = vector::zero;
            for (label pi = 0; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            // next find exact centre and the total solid angle, omega
            vector sumAc = vector::zero;
            scalar sumR = 0;
            scalar omega = 0;
            // The area vector is calculated as for Cartesian geometry
            // just to get the direction
            vector sumAn = vector::zero;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector c = p[f[pi]] + nextPoint + fCentre;
                scalar r = mag(p[f[pi]]) + mag(nextPoint);
                // The area vector as if for Cartesian geometry
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = sphTriSolidAngle(p[f[pi]], nextPoint, fCentre);
                omega += a;
                sumAc += a*c;
                sumR += r*a;
                sumAn += n;
            }
            if (mag(omega) < SMALL)
            {
                Info << "Face " << facei << " radii:\n";
                for(label pi = 0; pi < nPoints; pi++)
                {
                    Info << mag(p[f[pi]]) << ' ';
                }
                Info << "\nrhats:\n";
                for(label pi = 0; pi < nPoints; pi++)
                {
                    Info << p[f[pi]]/mag(p[f[pi]]) << ' ';
                }
                Info << endl;
            
                FatalErrorIn("makeFaceCentresAndAreas") << " for spherical geometry"
                    << " spherical face has a solid anlge of " << omega
                    << " and radii " << r1 << " and " << r2 << exit(FatalError);
            }
            
            sumR /= (2*omega);
            vector rhat = sumAc/mag(sumAc);
            fCtrs[facei] = rhat*sumR;
            magfAreas[facei] = sqr(sumR)*omega;
            fAreas[facei] = magfAreas[facei]*sign(sumAn & rhat)*rhat;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}


void Foam::primitiveMesh::overrideFaceCentres(const pointField& newCf)
{
    Info << "Overriding spherical face centres\n" << endl;
    
    if (!faceCentresPtr_)
    {
        calcFaceCentresAndAreas();
    }
    vectorField& faceCentres = *faceCentresPtr_;
    
    forAll(newCf, faceI)
    {
        faceCentres[faceI] = unitVector(newCf[faceI])*mag(faceCentres[faceI]);
    }
}


const Foam::vectorField& Foam::primitiveMesh::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::magFaceAreas() const
{
    if (!magFaceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *magFaceAreasPtr_;
}


// ************************************************************************* //
