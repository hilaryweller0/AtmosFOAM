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

#include "geostrophicZeroFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geostrophicZeroFluxFvPatchScalarField::
geostrophicZeroFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


geostrophicZeroFluxFvPatchScalarField::
geostrophicZeroFluxFvPatchScalarField
(
    const geostrophicZeroFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


geostrophicZeroFluxFvPatchScalarField::
geostrophicZeroFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary&
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


geostrophicZeroFluxFvPatchScalarField::
geostrophicZeroFluxFvPatchScalarField
(
    const geostrophicZeroFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf)
{}


geostrophicZeroFluxFvPatchScalarField::
geostrophicZeroFluxFvPatchScalarField
(
    const geostrophicZeroFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void geostrophicZeroFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& environmentalProperties
        = db().lookupObject<IOdictionary>("environmentalProperties");

    const dimensionedScalar g(environmentalProperties.lookup("g"));
    const dimensionedScalar beta(environmentalProperties.lookup("beta"));
    const dimensionedScalar f0(environmentalProperties.lookup("f0"));
    const vector OmegaHat(environmentalProperties.lookup("OmegaHat"));

    const surfaceVectorField& Uf
        = db().lookupObject<surfaceVectorField>("Uf");
    const fvsPatchField<vector> Ufp
        = patch().patchField<surfaceVectorField, vector>(Uf);

    const surfaceVectorField& Cf = Uf.mesh().Cf();
    fvsPatchField<vector> Cfp =
        patch().patchField<surfaceVectorField, vector>(Cf);

    gradient() = -(f0.value() + beta.value()*Cfp.component(vector::Y))*
                 ((OmegaHat ^ Ufp) & patch().nf())/g.value();
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void geostrophicZeroFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    geostrophicZeroFluxFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
