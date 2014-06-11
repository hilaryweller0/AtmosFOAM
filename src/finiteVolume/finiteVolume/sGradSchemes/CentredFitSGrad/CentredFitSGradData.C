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

#include "CentredFitSGradData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedCentredCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::CentredFitSGradData<Polynomial>::CentredFitSGradData
(
    const fvMesh& mesh,
    const extendedCentredCellToFaceStencil& stencil,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        CentredFitSGradData<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
    (
        mesh, stencil, true, linearLimitFactor, centralWeight
    ),
    delta_(mesh.delta()),
    coeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing CentredFitSGradData<Polynomial>" << endl;
    }

    calcFit();
    
//    Info << "CentredFitSGradData<Polynomial>::CentredFitSGradData() :"
//         << "coeffs = " << coeffs_ << endl;

    if (debug)
    {
        Info<< "CentredFitSGradData<Polynomial>::CentredFitSGradData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::CentredFitSGradData<Polynomial>::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const label facei
)
{
    const fvMesh& mesh = this->mesh();
    //const fvMeshWithDual& mesh
    //     = *(dynamic_cast<const fvMeshWithDual*>(&(this->mesh())));

    idir = delta_[facei];
    idir /= mag(idir);

    //if (mesh.gType() == fvMeshWithDual::CARTESIAN)
    // ONly for Cartesian geometry
    {
        if (mesh.nGeometricD() <= 2) // find the normal direction
        {
            if (mesh.geometricD()[0] == -1)
            {
                kdir = vector(1, 0, 0);
            }
            else if (mesh.geometricD()[1] == -1)
            {
                kdir = vector(0, 1, 0);
            }
            else
            {
                kdir = vector(0, 0, 1);
            }
        }
        else // 3D so find a direction in the plane of the face
        {
            const face& f = mesh.faces()[facei];
            kdir = mesh.points()[f[0]] - mesh.faceCentres()[facei];
        }
    }
//    else //Spherical geometry so kdir is the radial direction
//    {
//        kdir = mesh.faceCentres()[facei];
//    }

//    if (mesh.gType() != fvMeshWithDual::CARTESIAN || mesh.nGeometricD() == 3)
    if (mesh.nGeometricD() == 3)
    {
        // Remove the idir component from kdir and normalise
        kdir -= (idir & kdir)*idir;

        scalar magk = mag(kdir);

        if (magk < SMALL)
        {
            FatalErrorIn("findFaceDirs(..)") << " calculated kdir = zero"
                << exit(FatalError);
        }
        else
        {
            kdir /= magk;
        }
    }

    jdir = kdir ^ idir;
}


template<class Polynomial>
void Foam::CentredFitSGradData<Polynomial>::calcFit
(
    scalarList& coeffsi,
    const List<point>& C,
    const scalar wLin,
    const scalar deltaCoeff,
    const label facei
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, facei);

    // Setup the point weights
    scalarList wts(C.size(), scalar(1));
    wts[0] = this->centralWeight();
    wts[1] = this->centralWeight();

    // Reference point
    point p0 = this->mesh().faceCentres()[facei];

    // p0 -> p vector in the face-local coordinate system
    vector d;

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), this->minSize(), scalar(0));

    forAll(C, ip)
    {
        const point& p = C[ip];
        const vector p0p = p - p0;

        d.x() = p0p & idir;
        d.y() = p0p & jdir;
        d.z() = p0p & kdir;

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs(B[ip], d, wts[ip], this->dim());
    }

    // Additional weighting for constant and linear terms
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

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[1]*wts[i]*svd.VSinvUt()[1][i]/scale;
        }

        goodFit =
        (
            mag(wts[0]*wts[0]*svd.VSinvUt()[0][0] - wLin)
          < this->linearLimitFactor()*wLin)
         && (mag(wts[0]*wts[1]*svd.VSinvUt()[0][1] - (1 - wLin)
        ) < this->linearLimitFactor()*(1 - wLin))
         && coeffsi[0] < 0 && coeffsi[1] > 0
         && mag(coeffsi[0] + deltaCoeff) < 0.5*deltaCoeff
         && mag(coeffsi[1] - deltaCoeff) < 0.5*deltaCoeff;

        if (!goodFit)
        {
            // (not good fit so increase weight in the centre and weight
            //  for constant and linear terms)

            WarningIn
            (
                "CentredFitSGradData<Polynomial>::calcFit"
                "(const List<point>& C, const label facei"
            )   << "Cannot fit face " << facei << " iteration " << iIt
                << " with sum of weights " << sum(coeffsi) << nl
                << "    Weights " << coeffsi << nl
                << "    Linear weights " << wLin << " " << 1 - wLin << nl
                << "    deltaCoeff " << deltaCoeff << nl
                << "    sing vals " << svd.S() << nl
                << "Components of goodFit:\n"
                << "    wts[0]*wts[0]*svd.VSinvUt()[0][0] = "
                << wts[0]*wts[0]*svd.VSinvUt()[0][0] << nl
                << "    wts[0]*wts[1]*svd.VSinvUt()[0][1] = "
                << wts[0]*wts[1]*svd.VSinvUt()[0][1]
                << " dim = " << this->dim() << endl;

            wts[0] *= 10;
            wts[1] *= 10;

            for (label j = 0; j < B.m(); j++)
            {
                B[0][j] *= 10;
                B[1][j] *= 10;
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
        // Remove the uncorrected coefficients
        coeffsi[0] += deltaCoeff;
        coeffsi[1] -= deltaCoeff;
    }
    else
    {
        WarningIn
        (
            "CentredFitSGradData<Polynomial>::calcFit(..)"
        )   << "Could not fit face " << facei
            << "    Coefficients = " << coeffsi
            << ", reverting to uncorrected." << endl;

        coeffsi = 0;
    }
}


template<class Polynomial>
void Foam::CentredFitSGradData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles or tets.
    // Need bigger stencils
    List<List<point> > stencilPoints(mesh.nFaces());
    this->stencil().collectData(mesh.C(), stencilPoints);

    // find the fit coefficients for every face in the mesh

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const surfaceScalarField& dC = mesh.deltaCoeffs();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        calcFit
        (
            coeffs_[facei],
            stencilPoints[facei],
            w[facei],
            dC[facei],
            facei
        );
    }

    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();
    const surfaceScalarField::GeometricBoundaryField& bdC = dC.boundaryField();

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];
        const fvsPatchScalarField& pdC = bdC[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                calcFit
                (
                    coeffs_[facei],
                    stencilPoints[facei],
                    pw[i],
                    pdC[i],
                    facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //
