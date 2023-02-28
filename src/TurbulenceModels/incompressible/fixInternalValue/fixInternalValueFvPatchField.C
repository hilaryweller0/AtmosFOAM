/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fixInternalValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixInternalValueFvPatchField<Type>::fixInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    refValue_(p.size())
{}


template<class Type>
Foam::fixInternalValueFvPatchField<Type>::fixInternalValueFvPatchField
(
    const fixInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    zeroGradientFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(mapper((ptf.refValue_)))
{
    if (mappingRequired && notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::fixInternalValueFvPatchField<Type>::fixInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<Type>(p, iF, dict),
    refValue_("refValue", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::fixInternalValueFvPatchField<Type>::fixInternalValueFvPatchField
(
    const fixInternalValueFvPatchField& fivpf
)
:
    zeroGradientFvPatchField<Type>(fivpf),
    refValue_(fivpf.refValue())
{}


template<class Type>
Foam::fixInternalValueFvPatchField<Type>::fixInternalValueFvPatchField
(
    const fixInternalValueFvPatchField& fivpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(fivpf, iF),
    refValue_(fivpf.refValue())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixInternalValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    m(refValue_, refValue_);
}


template<class Type>
void Foam::fixInternalValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const fixInternalValueFvPatchField<Type>& mptf =
        refCast<const fixInternalValueFvPatchField<Type>>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
}


template<class Type>
void Foam::fixInternalValueFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    Field<Type>& internalField = const_cast<Field<Type>&>(this->primitiveField());

    Field<Type>::operator=(refValue_);
    
    forAll(*this, facei)
    {
        label celli = this->patch().faceCells()[facei];
        internalField[celli] = refValue_[facei];
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::fixInternalValueFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    // Apply the patch internal field as a constraint in the matrix
    matrix.setValues(this->patch().faceCells(), this->patchInternalField());
}


template<class Type>
void Foam::fixInternalValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "refValue", refValue_);
    writeEntry(os, "value", *this);
}

// ************************************************************************* //
