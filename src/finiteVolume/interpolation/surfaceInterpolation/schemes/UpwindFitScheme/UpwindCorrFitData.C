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
    const dictionary& debugDict = mesh.solutionDict().subOrEmptyDict("debug");

    // TODO: use this instead of debugFaceI
    labelList emptyList;
    labelList debugFaceIlist = debugDict.lookupOrDefault("interpolationWeightsFaceIndices", emptyList, false, false);

    label debugFaceI = -1;
    debugDict.readIfPresent("interpolationWeightsFaceIndex", debugFaceI, false, false);

    std::ostringstream ooss;
    ooss << "ownerWeights" << debugFaceI;

    volScalarField ownerWeights
    (
        IOobject
        (
            ooss.str(),
       	    mesh.time().constant(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
        ),
	mesh,
	0.0,
	"fixedValue"
    );

    std::ostringstream noss;
    noss << "neiWeights" << debugFaceI;

    volScalarField neiWeights
    (
        IOobject
        (
            noss.str(),
       	    mesh.time().constant(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
        ),
	mesh,
	0.0,
	"fixedValue"
    );


    // Get the cell/face centres in stencil order.
    List<List<point> > stencilPoints(mesh.nFaces());
    this->stencil().collectData
    (
        this->stencil().ownMap(),
        this->stencil().ownStencil(),
        mesh.C(),
        stencilPoints
    );

    fit(owncoeffs_, stencilPoints, ownerWeights, debugFaceI);
    ownerWeights.write();


    this->stencil().collectData
    (
        this->stencil().neiMap(),
        this->stencil().neiStencil(),
        mesh.C(),
        stencilPoints
    );

    fit(neicoeffs_, stencilPoints, neiWeights, debugFaceI);
    neiWeights.write();
}

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::fit(
        List<scalarList>& coeffs,
        const List<List<point> > stencilPoints,
        volScalarField& weightField,
        const label debugFaceI)
{
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        scalarList wts(stencilPoints[facei].size(), scalar(1));
        FitData
        <
            UpwindCorrFitData<Polynomial>,
            extendedUpwindCellToFaceStencilNew,
            Polynomial
        >::calcFit(coeffs[facei], wts, stencilPoints[facei], w[facei], facei);

        if (facei == debugFaceI)
        {
            forAll(stencilPoints[facei], stencilI)
            {
                forAll(mesh.C(), cellI)
                {
                    if (mesh.C()[cellI] == stencilPoints[facei][stencilI])
                    {
                        weightField[cellI] = coeffs[facei][stencilI];
                        if (stencilI == 0)
                        {
                            weightField[cellI] += 1.0; // re-add upwind weight
                        }
                    }
                }
                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& boundaryPatch = mesh.boundary()[patchI];
                    fvPatchField<scalar>& weightsPatch = weightField.boundaryField()[patchI];
                    forAll(boundaryPatch.Cf(), boundaryFaceI)
                    {
                        if (boundaryPatch.Cf()[boundaryFaceI] == stencilPoints[facei][stencilI])
                        {
                            weightsPatch[boundaryFaceI] = coeffs[facei][stencilI];
                            if (stencilI == 0)
                            {
                                weightsPatch[boundaryFaceI] += 1.0; // re-add upwind weight
                            }
                        }
                    }
                }
            }
        }
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
                scalarList wts(stencilPoints[facei].size(), scalar(1));
                FitData
                <
                    UpwindCorrFitData<Polynomial>,
                    extendedUpwindCellToFaceStencilNew,
                    Polynomial
                >::calcFit
                (
                    coeffs[facei], wts, stencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //
