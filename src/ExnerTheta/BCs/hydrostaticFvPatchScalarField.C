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

#include "hydrostaticFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrostaticFvPatchScalarField::
hydrostaticFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    gradPcoeffName_("gradPcoeff"),
    buoyancyName_("buoyancyf")
{}


hydrostaticFvPatchScalarField::
hydrostaticFvPatchScalarField
(
    const hydrostaticFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    gradPcoeffName_("gradPcoeff"),
    buoyancyName_("buoyancyf")
{}


hydrostaticFvPatchScalarField::
hydrostaticFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    gradPcoeffName_(dict.lookup("gradPcoeffName")),
    buoyancyName_(dict.lookup("buoyancyName"))
{
//    gradient() = 0;
//    fvPatchField<scalar>::operator=(patchInternalField());
}


hydrostaticFvPatchScalarField::
hydrostaticFvPatchScalarField
(
    const hydrostaticFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    gradPcoeffName_(wbppsf.gradPcoeffName_),
    buoyancyName_(wbppsf.buoyancyName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostaticFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& gradPcoeff =
        patch().lookupPatchField<surfaceScalarField, scalar>(gradPcoeffName_);

    const fvsPatchField<scalar>& buoyancyf =
        patch().lookupPatchField<surfaceScalarField, scalar>(buoyancyName_);

    gradient() = buoyancyf /(patch().magSf()*gradPcoeff);
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void hydrostaticFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
    Foam::writeEntry(os, "gradPcoeffName", gradPcoeffName_);
    Foam::writeEntry(os, "buoyancyName", buoyancyName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    hydrostaticFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
