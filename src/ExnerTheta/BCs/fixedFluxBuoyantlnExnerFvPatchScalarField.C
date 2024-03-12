/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "fixedFluxBuoyantlnExnerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedFluxBuoyantlnExnerFvPatchScalarField::
fixedFluxBuoyantlnExnerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


fixedFluxBuoyantlnExnerFvPatchScalarField::
fixedFluxBuoyantlnExnerFvPatchScalarField
(
    const fixedFluxBuoyantlnExnerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


fixedFluxBuoyantlnExnerFvPatchScalarField::
fixedFluxBuoyantlnExnerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{
//    gradient() = 0;
//    fvPatchField<scalar>::operator=(patchInternalField());
}


fixedFluxBuoyantlnExnerFvPatchScalarField::
fixedFluxBuoyantlnExnerFvPatchScalarField
(
    const fixedFluxBuoyantlnExnerFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedFluxBuoyantlnExnerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& environmentalProperties
        = db().lookupObject<IOdictionary>("environmentalProperties");

    const dictionary& thermoProperties
        = db().lookupObject<IOdictionary>("physicalProperties");

    dimensionedVector g(environmentalProperties.lookup("g"));

    const constTransport<hConstThermo<perfectGas<specie> > > air
    (
         "mixture", thermoProperties
    );

    const dimensionedScalar Cp("Cp", dimGasConstant, air.Cp(0,0));

    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    gradient() = (g.value() & patch().nf())/(Cp.value()*T);
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedFluxBuoyantlnExnerFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedFluxBuoyantlnExnerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
