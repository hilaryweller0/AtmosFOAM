/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fixedInternalValueFvPatchField

Description
    This boundary condition provides a mechanism to set boundary (cell) values
    directly into a matrix, i.e. to set a constraint condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixInternalValue;
        value           uniform 3.7;
    }
    \endverbatim

SourceFiles
    fixInternalValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixInternalValueFvPatchField_H
#define fixInternalValueFvPatchField_H

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fixInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixInternalValueFvPatchField
:
    public zeroGradientFvPatchField<Type>
{
    // Private Data

        //- Value field
        Field<Type> refValue_;


public:

    //- Runtime type information
    TypeName("fixInternalValue");


    // Constructors

        //- Construct from patch and internal field
        fixInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixInternalValueFvPatchField<Type>
        //  onto a new patch
        fixInternalValueFvPatchField
        (
            const fixInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Disallow copy without setting internal field reference
        fixInternalValueFvPatchField
        (
            const fixInternalValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixInternalValueFvPatchField
        (
            const fixInternalValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fixInternalValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //-Manipulate a matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
