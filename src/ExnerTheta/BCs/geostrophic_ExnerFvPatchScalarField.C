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

#include "geostrophic_ExnerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geostrophic_ExnerFvPatchScalarField::
geostrophic_ExnerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    Ug_(vector(0,0,0)),
    Omega_(vector(0,0,0))
{}


geostrophic_ExnerFvPatchScalarField::
geostrophic_ExnerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    Ug_(dict.lookup("geostrophicVelocity")),
    Omega_(dict.lookup("Omega"))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


geostrophic_ExnerFvPatchScalarField::
geostrophic_ExnerFvPatchScalarField
(
    const geostrophic_ExnerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    Ug_(ptf.Ug_),
    Omega_(ptf.Omega_)
{}


geostrophic_ExnerFvPatchScalarField::
geostrophic_ExnerFvPatchScalarField
(
    const geostrophic_ExnerFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    Ug_(wbppsf.Ug_),
    Omega_(wbppsf.Omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void geostrophic_ExnerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& thetaf
         = patch().lookupPatchField<surfaceScalarField, scalar>("thetaf");

    const dictionary& thermoProperties
        = db().lookupObject<IOdictionary>("physicalProperties");

    const constTransport<hConstThermo<perfectGas<specie> > > air
    (
         thermoProperties.subDict("mixture")
    );

    const dimensionedScalar Cp("Cp", dimGasConstant, air.Cp(0,0));
    gradient() = -2*((Omega_ ^ Ug_) & patch().nf())/(Cp.value()*thetaf);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void geostrophic_ExnerFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
    os << "        geostrophicVelocity " << Ug_ << ';' << endl;
    os << "        Omega               " << Omega_ << ';' << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    geostrophic_ExnerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
