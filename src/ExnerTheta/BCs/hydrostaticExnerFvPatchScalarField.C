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

#include "hydrostaticExnerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrostaticExnerFvPatchScalarField::
hydrostaticExnerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


hydrostaticExnerFvPatchScalarField::
hydrostaticExnerFvPatchScalarField
(
    const hydrostaticExnerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


hydrostaticExnerFvPatchScalarField::
hydrostaticExnerFvPatchScalarField
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


hydrostaticExnerFvPatchScalarField::
hydrostaticExnerFvPatchScalarField
(
    const hydrostaticExnerFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf)
{}


hydrostaticExnerFvPatchScalarField::
hydrostaticExnerFvPatchScalarField
(
    const hydrostaticExnerFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostaticExnerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& environmentalProperties
        = db().lookupObject<IOdictionary>("environmentalProperties");

    dimensionedVector g(environmentalProperties.lookup("g"));

    const fvsPatchField<scalar>& gradPcoeff =
        patch().lookupPatchField<surfaceScalarField, scalar>("gradPcoeff");

    gradient() = (g.value() & patch().nf())/gradPcoeff;
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void hydrostaticExnerFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    hydrostaticExnerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
