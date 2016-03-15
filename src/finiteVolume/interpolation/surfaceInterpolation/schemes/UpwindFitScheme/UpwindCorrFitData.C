/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "UpwindCorrFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "extendedUpwindCellToFaceStencilNew.H"
#include "StencilWeights.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::UpwindCorrFitData<Polynomial>::UpwindCorrFitData
(
    const fvMesh& mesh,
    const extendedUpwindCellToFaceStencilNew& stencil,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        UpwindCorrFitData<Polynomial>,
        extendedUpwindCellToFaceStencilNew,
        Polynomial
    >
    (
        mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
    ),
    owncoeffs_(mesh.nFaces()),
    neicoeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing UpwindCorrFitData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<< "UpwindCorrFitData<Polynomial>::UpwindCorrFitData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the cell/face centres in stencil order.
    List<List<point> > stencilPoints(mesh.nFaces());
    this->stencil().collectData
    (
        this->stencil().ownMap(),
        this->stencil().ownStencil(),
        mesh.C(),
        stencilPoints
    );

    StencilWeights ownerWeights(mesh, "ownerWeights");
    fit(owncoeffs_, stencilPoints, ownerWeights);
    ownerWeights.write();

    this->stencil().collectData
    (
        this->stencil().neiMap(),
        this->stencil().neiStencil(),
        mesh.C(),
        stencilPoints
    );

    StencilWeights neiWeights(mesh, "neiWeights");
    fit(neicoeffs_, stencilPoints, neiWeights);
    neiWeights.write();
}

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::fit(
        List<scalarList>& coeffs,
        const List<List<point> >& stencilPoints,
        StencilWeights& stencilWeights)
{
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        fit(facei, coeffs, stencilPoints, w[facei]);
        stencilWeights.fitted(facei, coeffs, stencilPoints);
    }

    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();
    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                fit(facei, coeffs, stencilPoints, pw[i]);
                facei++;
            }
        }
    }
}

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::fit(
        const label faceI,
        List<scalarList>& coeffs,
        const List<List<point> >& stencilPoints,
        const scalar wLin)
{
    scalarList wts(stencilPoints[faceI].size(), scalar(1));
    autoPtr<Fit> fit = FitData
    <
        UpwindCorrFitData<Polynomial>,
        extendedUpwindCellToFaceStencilNew,
        Polynomial
    >::calcFit
    (
        coeffs[faceI], wts, stencilPoints[faceI], wLin, faceI
    );
    if (!fit->good)
    {
        WarningIn
        (
            "FitData<Polynomial>::calcFit(..)"
        )   << "Could not fit face " << faceI
            << "    Weights = " << coeffs[faceI]
            << ", reverting to upwind/linear." << nl
            << "    Linear weights " << wLin << " " << 1 - wLin << endl;
    }
}
