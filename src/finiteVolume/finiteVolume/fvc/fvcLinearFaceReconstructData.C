/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fvcLinearFaceReconstructData.H"
#include "surfaceFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "centredCFCFaceToFaceStencilObject.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fvcLinearFaceReconstructData::fvcLinearFaceReconstructData
(
    const fvMesh& mesh,
    const centredCFCFaceToFaceStencilObject& stencil
)
:
    MeshObject<fvMesh, fvcLinearFaceReconstructData>(mesh),
    stencil_(stencil),
    dim_(mesh.nGeometricD()),
    spherical_
    (
        mesh.nGeometricD() != cmptSum(mesh.geometricD() + Vector<label>::one)/2
    ),
    nValidFaces_(mesh.nFaces() - (3-mesh.nGeometricD())*2*mesh.nCells()),
    coeffs_(nValidFaces_)
{
    calcFit();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvcLinearFaceReconstructData::addCoeffs
(
    scalar* coeffs,
    const vector& d,
    const vector& Sf,
    const label faceI,
    const bool smallStencil
)
{
    // Matrix coefficients for flux terms
    //const fvMesh& mesh = this->mesh();
    
    const scalar w = magSqr(d) < SMALL ? 1e3 : 1.0;
    const scalar w2 = 1e3;

    label ic = 0;
    coeffs[ic++] = w*w2*Sf.x();                          // ax
    coeffs[ic++] = w*Sf.y();                          // ay
    if (dim_ == 3) coeffs[ic++] = w*Sf.z();           // az

    coeffs[ic++] = w*d.x()*Sf.x();                    // bx
    
    if (!smallStencil)
    {
        coeffs[ic++] = w*d.x()*Sf.y();                // by
        if (dim_ == 3) coeffs[ic++] = w*d.x()*Sf.z(); // bz

        //coeffs[ic++] = d.y()*Sf.x();                // cx
        coeffs[ic++] = w*d.y()*Sf.y();                // cy
        
        if (dim_ == 3)
        {
            //coeffs[ic++] = d.y()*Sf.z();      // cz
            //coeffs[ic++] = d.z()*Sf.x();      // dx
            //coeffs[ic++] = d.z()*Sf.y();      // dy
            coeffs[ic++] = w*d.z()*Sf.z();      // dz
        }
    }
}

void Foam::fvcLinearFaceReconstructData::findFaceDirs
(
    vector& idir,         // value changed in return
    vector& jdir,         // value changed in return
    vector& kdir,         // value changed in return
    const label faceI
)
{
    const fvMesh& mesh = this->mesh();
    const pointField& Cf = mesh.faceCentres();
    const pointField& C = mesh.cellCentres();
    const point C0 = Cf[faceI];
    const labelList& cellCells = mesh.cellCells()[mesh.faceOwner()[faceI]];

    // Set local i, j and k directions
    // idir is normal to the face
    idir = mesh.faceAreas()[faceI];
    idir /= mag(idir);
    
    // For spherical geometry kdir is the radial direction
    if (spherical_)
    {
        kdir = Cf[faceI];
        kdir -= (kdir & idir)*idir;
        kdir /= mag(kdir);
        jdir = kdir ^ idir;
    }
    else
    {
        // jdir is a direction of the cellCells of the owner of the face
        scalar minDotProd = 1;
        for(label i = 0; i < cellCells.size(); i++)
        {
            vector d = C[cellCells[i]] - C0;
            d /= mag(d);
            if (mag(idir & d) < minDotProd)
            {
                jdir = d;
                minDotProd = mag(idir & d);
            }
        }
        // kdir is out of the face of idir and jdir
        kdir = idir ^ jdir;
        kdir /= mag(kdir);
        // make sure jdir is at right angles to idir
        jdir -= (idir & jdir) * idir;
        jdir /= mag(jdir);
    }
    
    // Check orthogonality and size of coord system
    if
    (
        mag(idir & jdir) > SMALL
     || mag(idir & kdir) > SMALL
     || mag(jdir & kdir) > SMALL
     || mag(mag(idir) - 1) > SMALL
     || mag(mag(jdir) - 1) > SMALL
     || mag(mag(kdir) - 1) > SMALL
    )
    {
        FatalErrorIn("Foam::fvcLinearFaceReconstructData::findFaceDirs")
            << " face " << faceI << " has incorrect coord system with:\n"
            << "idir = " << idir << nl
            << "jdir = " << jdir << nl
            << "kdir = " << kdir << nl
            << "mag(idir & jdir) = " << mag(idir & jdir) << nl
            << "mag(idir & kdir) = " << mag(idir & kdir) << nl
            << "mag(jdir & jkir) = " << mag(jdir & kdir) << nl
            << "mag(idir) - 1 = " << mag(idir) - 1 << nl
            << "mag(jdir) - 1 = " << mag(jdir) - 1 << nl
            << "mag(kdir) - 1 = " << mag(kdir) - 1 << endl
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvcLinearFaceReconstructData::calcFit
(
    List<vector>& coeffsi,
    const List<point>& Cf,
    const List<vector>& Sf,
    const label faceI
)
{
    const fvMesh& mesh = this->mesh();

    const label stencilSize = Cf.size();
    coeffsi.setSize(stencilSize);

    // Local coord system for this face
    vector idir, jdir, kdir;

    findFaceDirs(idir, jdir, kdir, faceI);

    tensor coordTransform
    (
        idir.x(), idir.y(), idir.z(),
        jdir.x(), jdir.y(), jdir.z(),
        kdir.x(), kdir.y(), kdir.z()
    );
    tensor invTransform = inv(coordTransform);

//   Polynomial fit
    // Central face point
    point p0 = mesh.faceCentres()[faceI];
    //point p0 = this->stencil().origin(faceI);

    // Determine if the possible number of terms to be fit from the size
    // of the stencil
    const bool smallStencil = dim_ == 2 ?
        stencilSize <= 5 ? true : false
      : stencilSize <= 10 ? true : false;

    // Data local to each face:
    const label minSize = dim_ == 2 ?
          smallStencil ? 3 : 5
        : smallStencil ? 4 : 8;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(Cf.size(), minSize, scalar(0));

    for(label ip = 0; ip < Cf.size(); ip++)
    {
        // create vectors d and S in the local co-ordinate system
        vector d = Cf[ip] - p0;
        vector S = Sf[ip];
        
//        if (!spherical_)
//        {
            d = vector(d&idir, d&jdir, d&kdir);
            S = vector(S&idir, S&jdir, S&kdir);
//        }
//        else
//        {
//            // Local coord system at midpoint of d
//            point midd = 0.5*(Cf[ip] + p0);
//            vector k(midd/mag(midd));
//            vector i(idir - (idir & k)*k);
//            i /= mag(i);
//            vector j(k^i);

//            d = vector(d&i, d&j, d&k);

//            // Local coord system at Sf
//            k = Cf[ip]/mag(Cf[ip]);
//            i = idir - (idir & k)*k;
//            i /= mag(i);
//            j = k^i;

//            S = vector(S&i, S&j, S&k);
//        }

        // Add the polynomial coefficients to the matrix
        addCoeffs(B[ip], d, S, faceI, smallStencil);
    }

    // Calcultate the fit
    SVD svd(B, SMALL);

    vector a(1/(Sf[0]&idir),0,0);

    for(label i=1; i<stencilSize; i++)
    {
        a[1] -= a[0]*svd.VSinvUt()[1][i]*(Sf[i]&idir);
        if (dim_ == 3) a[2] -= a[0]*svd.VSinvUt()[2][i]*(Sf[i]&idir);
    }
    
    coeffsi[0] = invTransform & a;

    for(label i=1; i<stencilSize; i++)
    {
        vector a(0,0,0);

        a[1] = svd.VSinvUt()[1][i];
        if (dim_ == 3) a[2] = svd.VSinvUt()[2][i];

        coeffsi[i] = invTransform & a;
    }
}

void Foam::fvcLinearFaceReconstructData::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the face centres and area vectors in stencil order.
    List<List<point> > stencilPoints(nValidFaces_);
    this->stencil().collectData(mesh.Cf(), stencilPoints);
    List<List<vector> > stencilSf(nValidFaces_);
    this->stencil().collectData(mesh.Sf(), stencilSf);

    // find the fit coefficients for every internal face in the mesh

    for(label faceI = 0; faceI < nValidFaces_; faceI++)
    {
        fvcLinearFaceReconstructData::calcFit
        (
            coeffs_[faceI],
            stencilPoints[faceI],
            stencilSf[faceI],
            faceI
        );
    }
}


bool Foam::fvcLinearFaceReconstructData::movePoints()
{
    calcFit();
    return true;
}

// ************************************************************************* //
