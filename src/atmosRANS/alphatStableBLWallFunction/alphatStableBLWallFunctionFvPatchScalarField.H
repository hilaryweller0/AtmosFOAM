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

Class
    Foam::alphatStableBLWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent thermal diffusivity condition
    based on stableBLkEpsilon.

    Reference:
    \verbatim
        Cuxart et al, (2006) Boundary Layer Meteorology
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            alphatStableBLWallFunction;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::alphatWallFunctionFvPatchScalarField

SourceFiles
    alphatStableBLWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatStableBLWallFunctionFvPatchScalarField_H
#define alphatStableBLWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class alphatStableBLWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatStableBLWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("alphatStableBLWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatStableBLWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatStableBLWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatStableBLWallFunctionFvPatchScalarField
        (
            const alphatStableBLWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatStableBLWallFunctionFvPatchScalarField
        (
            const alphatStableBLWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatStableBLWallFunctionFvPatchScalarField
        (
            const alphatStableBLWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatStableBLWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
