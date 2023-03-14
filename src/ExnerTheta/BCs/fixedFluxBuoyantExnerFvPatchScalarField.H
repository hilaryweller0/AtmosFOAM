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

Class
    Foam::fixedFluxBuoyantExnerFvPatchScalarField

Description
    Foam::fixedFluxBuoyantExnerFvPatchScalarField

SourceFiles
    fixedFluxBuoyantExnerFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedFluxBuoyantExnerFvPatchScalarFields_H
#define fixedFluxBuoyantExnerFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedFluxBuoyantExnerFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedFluxBuoyantExnerFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("fixedFluxBuoyantExner");


    // Constructors

        //- Construct from patch and internal field
        fixedFluxBuoyantExnerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedFluxBuoyantExnerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedFluxBuoyantExnerFvPatchScalarField onto a new patch
        fixedFluxBuoyantExnerFvPatchScalarField
        (
            const fixedFluxBuoyantExnerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedFluxBuoyantExnerFvPatchScalarField
        (
            const fixedFluxBuoyantExnerFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        fixedFluxBuoyantExnerFvPatchScalarField
        (
            const fixedFluxBuoyantExnerFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFluxBuoyantExnerFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
