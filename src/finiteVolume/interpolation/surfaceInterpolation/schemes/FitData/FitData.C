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
#include "PolynomialFit.H"
#include "FixedPolynomial.H"
#include "Basis.H"
#include "surfaceFields.H"
#include "volFields.H"

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
    dim_
    (
        mesh.nGeometricD() == 1 ? 1 :
        isA<emptyPolyPatch>(mesh.boundaryMesh().last()) ? 2 :
        mesh.nGeometricD()
    )
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

    PolynomialFit<FixedPolynomial<Polynomial> > polynomialFit(
        linearCorrection_,
        linearLimitFactor_,
        centralWeight_, 
        dim_
    );
    const Basis basis(idir, jdir, kdir);
    polynomialFit.fit(coeffsi, wts, C, wLin, p0, pureUpwind, basis, facei);
}

template<class FitDataType, class ExtendedStencil, class Polynomial>
bool Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::movePoints()
{
    calcFit();
    return true;
}

// ************************************************************************* //
