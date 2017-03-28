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
#include "Basis.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fitCoefficients.H"
#include "localStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Form, class ExtendedStencil, class Polynomial>
Foam::FitData<Form, ExtendedStencil, Polynomial>::FitData
(
    const fvMesh& mesh,
    const ExtendedStencil& stencil,
    const bool linearCorrection,
    const scalar centralWeight,
    const bool sphericalGeometry
)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, Form>(mesh),
    stencil_(stencil),
    linearCorrection_(linearCorrection),
    centralWeight_(centralWeight),
    dim_
    (
        mesh.nGeometricD() == 1 ? 1 :
        mesh.nGeometricD() == 0 ? 2 :
        isA<emptyPolyPatch>(mesh.boundaryMesh().last()) ? 2 :
        mesh.nGeometricD()
    ),
    sphericalGeometry_(sphericalGeometry)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FitDataType, class ExtendedStencil, class Polynomial>
void Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const label facei,
    const List<point>& C
)
{
    const fvMesh& mesh = this->mesh();

    idir = mesh.faceAreas()[facei];
    idir /= mag(idir);
    
    const point& fC = mesh.faceCentres()[facei];

    if (sphericalGeometry_)
    {
        kdir = fC/mag(fC);
        jdir = kdir ^ idir;
    }
    else if (dim_ == 3)
    {
        List<vector> unitVectors(3);
        unitVectors[0] = vector(0, 0, 1);
        unitVectors[1] = vector(0, 1, 0);
        unitVectors[2] = vector(1, 0, 0);

        jdir = idir;
        forAll(unitVectors, vi)
        {
            if (magSqr(unitVectors[vi] ^ idir) > magSqr(jdir ^ idir))
            {
                jdir = unitVectors[vi];
            }
        }
        
        // Remove the idir from jdir and then normalise jdir
        jdir = jdir - (idir & jdir) * idir;
        jdir /= mag(jdir);
        
        // kdir is normal to idir and jdir
        kdir = idir ^ jdir;
    }
    else
    {
        // For the jdir find the cell within the stencil that gives the
        // direction most different to idir by finding the smallest idir&jdir
        jdir = idir;
        vector jdirTmp;

        forAll(C, celli)
        {
            jdirTmp = (C[celli] - fC)/mag(C[celli] - fC);
            if (mag(jdirTmp & idir) < mag(jdir & idir))
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
}

template<class FitDataType, class ExtendedStencil, class Polynomial>
autoPtr<fitResult> Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::calcFit
(
    fitCoefficients& coefficients,
    const List<point>& C,
    const label facei,
    bool owner
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, facei, C);
    
    point p0 = this->mesh().faceCentres()[facei];

    scalar c0SurfaceNormalComponent = this->mesh().faceAreas()[facei]
                                    & (C[0]-p0);
    scalar c1SurfaceNormalComponent = this->mesh().faceAreas()[facei]
                                    & (C[1]-p0);
    bool pureUpwind = 
    (
        sign(c0SurfaceNormalComponent) == sign(c1SurfaceNormalComponent)
    );

    fitWeights weights(C.size());
    weights.setCentralWeight(centralWeight_, pureUpwind);

    PolynomialFit<Polynomial> polynomialFit(dim_, 1e-9);

    const Basis basis(idir, jdir, kdir);
    const localStencil stencil(C, p0, basis);

    return polynomialFit.fit(coefficients, weights, stencil);
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
    fitCoefficients coefficients(C.size(), linearCorrection_, wLin);
    calcFit(coefficients, C, facei);
    coefficients.copyInto(coeffsi);
    // FIXME: I should probably populate wts using some returned data inside fitResult
}

template<class FitDataType, class ExtendedStencil, class Polynomial>
bool Foam::FitData<FitDataType, ExtendedStencil, Polynomial>::movePoints()
{
    calcFit();
    return true;
}

// ************************************************************************* //
