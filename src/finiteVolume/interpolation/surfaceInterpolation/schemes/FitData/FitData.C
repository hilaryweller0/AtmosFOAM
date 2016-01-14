/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "FitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Form, class ExtendedStencil, class Polynomial>
Foam::FitData<Form, ExtendedStencil, Polynomial>::FitData
(
    const fvMesh& mesh,
    const ExtendedStencil& stencil,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, Form>(mesh),
    stencil_(stencil),
    linearCorrection_(linearCorrection),
    linearLimitFactor_(linearLimitFactor),
    centralWeight_(centralWeight),
// 2D dependent on the size of the 
    dim_
    (
        mesh.nGeometricD() == 1 ? 1 :
        isA<emptyPolyPatch>(mesh.boundaryMesh().last()) ? 2 :
        mesh.nGeometricD()
    ),
    minSize_(Polynomial::nTerms(dim_))
{
    // Check input
    if (linearLimitFactor <= SMALL || linearLimitFactor > 4)
    {
        FatalErrorIn("FitData<Polynomial>::FitData(..)")
            << "linearLimitFactor requested = " << linearLimitFactor
            << " should be between zero and 4"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FitDataType, class ExtendedStencil, class Polynomial>
void Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const label facei
)
{
    const fvMesh& mesh = this->mesh();

    idir = mesh.faceAreas()[facei];
    idir /= mag(idir);
    
    const point& fC = mesh.faceCentres()[facei];
    
    // For the jdir find the cell within all cellCells of the owner that gives
    // the direction most different to idir by finding the largest idir ^ jdir
    jdir = idir;
    vector jdirTmp;
    const label own = mesh.faceOwner()[facei];
    forAll(mesh.cellCells()[own], i)
    {
        const label celli = mesh.cellCells()[own][i];
        jdirTmp = mesh.cellCentres()[celli] - fC;
        if (magSqr(jdirTmp ^ idir) > magSqr(jdir ^ idir))
        {
            jdir = jdirTmp;
        }
    }
    // Remove the idir from jdir and then normalise jdir
    jdir = jdir - (idir & jdir) * idir;
    jdir /= mag(jdir);
    
    // kdir is normal to idir and jdir
    kdir = idir ^ jdir;
}


template<class FitDataType, class ExtendedStencil, class Polynomial>
void Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::calcFit
(
    scalarList& coeffsi,
    scalarList& wts,
    const List<point>& C,
    const scalar wLin,
    const label facei
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, facei);
    
    // Reference point
    point p0 = this->mesh().faceCentres()[facei];

    scalar c0SurfaceNormalComponent = this->mesh().faceAreas()[facei] & (C[0]-p0);
    scalar c1SurfaceNormalComponent = this->mesh().faceAreas()[facei] & (C[1]-p0);
    bool pureUpwind = (sign(c0SurfaceNormalComponent) == sign(c1SurfaceNormalComponent));

    // Setup the point weights
    wts[0] = centralWeight_;
    if (!pureUpwind)
    {
        wts[1] = centralWeight_;
    }

    // p0 -> p vector in the face-local coordinate system
    vector d;

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), minSize_, scalar(0));

    forAll(C, ip)
    {
        const point& p = C[ip];

        d.x() = (p - p0)&idir;
        d.y() = (p - p0)&jdir;
        #ifndef SPHERICAL_GEOMETRY
        d.z() = (p - p0)&kdir;
        #else
        d.z() = mag(p) - mag(p0);
        #endif

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs(B[ip], d, wts[ip], dim_);
    }

    // Additional weighting for constant (and linear) terms
    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= wts[0];
        B[i][1] *= wts[0];
    }

    // Set the fit
    label stencilSize = C.size();
    coeffsi.setSize(stencilSize);

    bool goodFit = false;
    for (int iIt = 0; iIt < 8 && !goodFit; iIt++)
    {
        SVD svd(B, SMALL);

        scalar maxCoeff = 0;
        label maxCoeffi = 0;

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[0]*wts[i]*svd.VSinvUt()[0][i];
            if (mag(coeffsi[i]) > maxCoeff)
            {
                maxCoeff = mag(coeffsi[i]);
                maxCoeffi = i;
            }
        }

        if (linearCorrection_)
        {
            goodFit =
                (mag(coeffsi[0] - wLin) < linearLimitFactor_*wLin)
             && (mag(coeffsi[1] - (1 - wLin)) < linearLimitFactor_*(1 - wLin))
             && maxCoeffi <= 1;
        }
        else
        {
            // go through all coeffsi except the first, add up all positive coeffs
            scalar positiveCoeffSum = 0;
            forAll(coeffsi, i)
            {
                if (i > 0 && coeffsi[i] > 0)
                {
                    positiveCoeffSum += coeffsi[i];
                }
            }

            // Upwind: weight on face is 0
            goodFit = (mag(coeffsi[0] - 1.0) < linearLimitFactor_*1.0)
                && coeffsi[0] > positiveCoeffSum;
        }

        if (!goodFit) // (not good fit so increase weight in the centre and
                      //  weight for constant and linear terms)
        {

            wts[0] *= 10;
            if (linearCorrection_)
            {
                wts[1] *= 10;
            }
            else if (!pureUpwind && iIt == 0)
            {
                wts[1] /= centralWeight_;
            }

            for (label j = 0; j < B.m(); j++)
            {
                B[0][j] *= 10;
                if (linearCorrection_) B[1][j] *= 10;
                else if (!pureUpwind && iIt == 0) B[1][j] /= centralWeight_;
            }

            for (label i = 0; i < B.n(); i++)
            {
                B[i][0] *= 10;
                B[i][1] *= 10;
            }
        }
    }

    if (goodFit)
    {
        if (linearCorrection_)
        {
            // Remove the uncorrected linear coefficients
            coeffsi[0] -= wLin;
            coeffsi[1] -= 1 - wLin;
        }
        else
        {
            // Remove the uncorrected upwind coefficients
            coeffsi[0] -= 1.0;
        }
    }
    else
    {
        WarningIn
        (
            "FitData<Polynomial>::calcFit(..)"
        )   << "Could not fit face " << facei
            << "    Weights = " << coeffsi
            << ", reverting to upwind/linear." << nl
            << "    Linear weights " << wLin << " " << 1 - wLin << endl;

        coeffsi = 0;
        
        if (linearCorrection_)
        {
            coeffsi[0] = 1-wLin;
            coeffsi[1] = -(1-wLin);
        }
    }
}


template<class FitDataType, class ExtendedStencil, class Polynomial>
bool Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::movePoints()
{
    calcFit();
    return true;
}

// ************************************************************************* //
