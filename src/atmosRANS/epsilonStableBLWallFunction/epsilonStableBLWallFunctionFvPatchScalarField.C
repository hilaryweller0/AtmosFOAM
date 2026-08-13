/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "epsilonStableBLWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::epsilonStableBLWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& mtm
) const
{
    const label patchi = patch().index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(mtm, patchi);

    const scalarField& y = mtm.yb()[patchi];

    const tmp<scalarField> tnuw = mtm.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = mtm.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = mtm.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const volScalarField& epsilonGroundCorr
         = db().lookupObject<volScalarField>("epsilonGroundCorr");
    const scalarField& groundCorr = epsilonGroundCorr.boundaryField()[patchi];

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu75 = pow(nutw.Cmu(), 0.75);

    // Initialise the returned field pair
    Pair<tmp<scalarField>> GandEpsilon
    (
        tmp<scalarField>(new scalarField(patch().size())),
        tmp<scalarField>(new scalarField(patch().size()))
    );
    scalarField& G = GandEpsilon.first().ref();
    scalarField& epsilon  = GandEpsilon.second().ref();

    // Calculate G and epsilon
    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];

        G[facei] =
            (nutw[facei] + nuw[facei])
           *magGradUw[facei]
           *Cmu25*sqrt(k[celli])
          /(nutw.kappa()*y[facei]);

        epsilon[facei] = Cmu75*k[celli]*sqrt(k[celli])/(nutw.kappa()*y[facei])
                        *groundCorr[facei];
        Info << "epsilon ground corr = " << groundCorr[facei] << nl;
    }

    return GandEpsilon;
}


void Foam::epsilonStableBLWallFunctionFvPatchScalarField::updateCoeffsMaster()
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    // Lookup the momentum transport model
    const momentumTransportModel& mtm =
        db().lookupType<momentumTransportModel>(internalField().group());

    // Get mutable references to the turbulence fields
    VolInternalField<scalar>& G =
        db().lookupObjectRef<VolInternalField<scalar>>(mtm.GName());
    volScalarField& epsilon =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    const volScalarField::Boundary& epsilonBf = epsilon.boundaryField();

    // Make all processors build the near wall distances
    mtm.yb();

    // Evaluate all the wall functions
    PtrList<scalarField> Gpfs(epsilonBf.size());
    PtrList<scalarField> epsilonPfs(epsilonBf.size());
    forAll(epsilonBf, patchi)
    {
        if (isA<epsilonStableBLWallFunctionFvPatchScalarField>(epsilonBf[patchi]))
        {
            const epsilonStableBLWallFunctionFvPatchScalarField& epsilonPf =
                refCast<const epsilonStableBLWallFunctionFvPatchScalarField>
                (epsilonBf[patchi]);

            Pair<tmp<scalarField>> GandEpsilon = epsilonPf.calculate(mtm);

            Gpfs.set(patchi, GandEpsilon.first().ptr());
            epsilonPfs.set(patchi, GandEpsilon.second().ptr());
        }
    }

    // Average the values into the wall-adjacent cells and store
    wallCellGPtr_.reset(patchFieldsToWallCellField(Gpfs).ptr());
    wallCellEpsilonPtr_.reset(patchFieldsToWallCellField(epsilonPfs).ptr());

    // Set the fractional proportion of the value in the wall cells
    UIndirectList<scalar>(G, wallCells()) =
        (1 - wallCellFraction())*scalarField(G, wallCells())
      + wallCellFraction()*wallCellGPtr_();
    UIndirectList<scalar>(epsilon, wallCells()) =
        (1 - wallCellFraction())*scalarField(epsilon, wallCells())
      + wallCellFraction()*wallCellEpsilonPtr_();
}


void Foam::epsilonStableBLWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues(wallCells(), wallCellEpsilonPtr_(), wallCellFraction());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonStableBLWallFunctionFvPatchScalarField::
epsilonStableBLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


Foam::epsilonStableBLWallFunctionFvPatchScalarField::
epsilonStableBLWallFunctionFvPatchScalarField
(
    const epsilonStableBLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


Foam::epsilonStableBLWallFunctionFvPatchScalarField::
epsilonStableBLWallFunctionFvPatchScalarField
(
    const epsilonStableBLWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(ewfpsf, iF),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::epsilonStableBLWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    wallCellWallFunctionFvPatchScalarField::map(ptf, mapper);
    wallCellGPtr_.clear();
    wallCellEpsilonPtr_.clear();
}


void Foam::epsilonStableBLWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    wallCellWallFunctionFvPatchScalarField::reset(ptf);
    wallCellGPtr_.clear();
    wallCellEpsilonPtr_.clear();
}


void Foam::epsilonStableBLWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    initMaster();

    updateCoeffsMaster();

    operator==(patchInternalField());

    fvPatchScalarField::updateCoeffs();
}


void Foam::epsilonStableBLWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    if (masterPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "updateCoeffs must be called before manipulateMatrix"
            << exit(FatalError);
    }

    manipulateMatrixMaster(matrix);

    fvPatchScalarField::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonStableBLWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
