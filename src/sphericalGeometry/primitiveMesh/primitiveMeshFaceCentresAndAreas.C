/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

    For spherical geometry - only copes with either purely radial or purely
    spherical faces - not a mixture. A spherical face is a face on the sphere
    whereas a radial face points around the sphere

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const scalar RS_TOL = 1e-7;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcFaceCentresAndAreas() const
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
    if (faceCentresPtr_ || faceAreasPtr_)
    {
        FatalErrorIn("primitiveMesh::calcFaceCentresAndAreas() const")
            << "Face centres or face areas already calculated"
            << exit(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas);

    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Finished calculating face centres and face areas"
            << endl;
    }
}


void primitiveMesh::makeFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas
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
                const vector rhat = 0.5*(rhat0 + rhat1);
                fCtrs[facei] = 0.5*(r1 + r2)*rhat;
                const vector idir = rhat ^ (rhat1 - rhat0);
                //const scalar A = arcLength(rhat0, rhat1)*Rsphere*(r2 - r1);
                const scalar A = 0.5*arcLength(rhat0, rhat1)*mag(r2S - r1S);
                fAreas[facei] = A*idir/mag(idir);
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
            vector sumAn = vector::zero;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector c = p[f[pi]] + nextPoint + fCentre;
                scalar r = mag(p[f[pi]]) + mag(nextPoint);
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
            //fAreas[facei] = RsphereS*rhat*omega;
            fAreas[facei] = sqr(sumR)*rhat*omega*sign(sumAn & rhat);
            
        }
    }

//    // Test that face area vectors of each cell sum to zero
//    forAll(cells(), cellI)
//    {
//        const cell& c = cells()[cellI];
//        vector sumS = vector::zero;
//        scalar maxFaceAreaS = 0;
//        forAll(c, facei)
//        {
//            if (cellI == faceOwner()[c[facei]]) sumS += fAreas[c[facei]];
//            else sumS -= fAreas[c[facei]];
//            
//            scalar faceAreaS = magSqr(fAreas[c[facei]]);
//            if (faceAreaS > maxFaceAreaS)
//            {
//                maxFaceAreaS = faceAreaS;
//            }
//        }
//    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& primitiveMesh::faceCentres() const
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


const vectorField& primitiveMesh::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
