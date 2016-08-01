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
#include "stencilWeights.H"

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
        Info<< "Constructing UpwindCorrFitData<Polynomial>" << endl;
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

    Info << "===OWNER===" << endl;
    stencilWeights ownerWeights(mesh, "owner");
    fit(owncoeffs_, stencilPoints, ownerWeights, true);
    ownerWeights.write();

    this->stencil().collectData
    (
        this->stencil().neiMap(),
        this->stencil().neiStencil(),
        mesh.C(),
        stencilPoints
    );

    Info << "===NEIGHBOUR===" << endl;
    stencilWeights neiWeights(mesh, "nei");
    fit(neicoeffs_, stencilPoints, neiWeights, false);
    neiWeights.write();
}

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::fit(
        List<scalarList>& coeffs,
        const List<List<point> >& stencilPoints,
        stencilWeights& stencilWts,
        bool owner
)
{
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        fitCoefficients coefficients
        (
            stencilPoints[facei].size(),
            FitData<
                UpwindCorrFitData<Polynomial>,
                extendedUpwindCellToFaceStencilNew,
                Polynomial
            >::linearCorrection(),
            w[facei]
        );
        autoPtr<fitResult> f = fit(facei, coefficients, stencilPoints[facei], owner);
        coefficients.copyInto(coeffs[facei]);
        stencilWts.fitted(facei, f, stencilPoints[facei]);
    }

    const surfaceScalarField::Boundary& bw = w.boundaryField();
    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                fitCoefficients coefficients
                (
                    stencilPoints[facei].size(),
                    FitData<
                        UpwindCorrFitData<Polynomial>,
                        extendedUpwindCellToFaceStencilNew,
                        Polynomial
                    >::linearCorrection(),
                    pw[i]
                );
                autoPtr<fitResult> f = fit(facei, coefficients, stencilPoints[facei], owner);
                coefficients.copyInto(coeffs[facei]);
                stencilWts.fitted(facei, f, stencilPoints[facei]);
                facei++;
            }
        }
    }
}

template<class Polynomial>
autoPtr<fitResult> Foam::UpwindCorrFitData<Polynomial>::fit
(
    const label faceI,
    fitCoefficients& coefficients,
    const List<point>& stencilPoints,
    bool owner
)
{
    autoPtr<fitResult> fit = FitData
    <
        UpwindCorrFitData<Polynomial>,
        extendedUpwindCellToFaceStencilNew,
        Polynomial
    >::calcFit
    (
        coefficients, stencilPoints, faceI, owner
    );

    if (!fit->good)
    {
/*        WarningIn
        (
            "FitData<Polynomial>::calcFit(..)"
        )   << "Could not fit face " << faceI
            << " at " << this->mesh().Cf()[faceI]
            << "    Weights = " << coefficients
            << ", reverting to upwind/linear." << endl;
            */
    }

    return fit;
}
