/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "PartitionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const IOobject& io,
    const Mesh& mesh,
    const PartitionedFraction<Type, PatchField, GeoMesh>& sigma__
)
:
    PartitionedFraction<Type, PatchField, GeoMesh>(partNames__, io, mesh),
    sigma_(sigma__)
{
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const PartitionedFraction<Type, PatchField, GeoMesh>& sigma__
)
:
    PartitionedFraction<Type, PatchField, GeoMesh>(partNames__, field),
    sigma_(sigma__)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const PartitionedField<Type, PatchField, GeoMesh>& copyFrom
)
:
    PartitionedFraction<Type, PatchField, GeoMesh>(copyFrom),
    sigma_(copyFrom.sigma())
{}

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& newName,
    const PartitionedField<Type, PatchField, GeoMesh>& copyFrom
)
:
    PartitionedFraction<Type, PatchField, GeoMesh>(copyFrom, newName),
    sigma_(copyFrom.sigma())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::~PartitionedField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::PartitionedField<Type, PatchField, GeoMesh>::updateSum()
{
    GeometricField<Type, PatchField, GeoMesh>& sum_
         = PartitionedFraction<Type, PatchField, GeoMesh>::sum_;

    sum_ = sigma_[0]*operator[](0);

    // Sum contributions from other partitions
    for (label ip = 1; ip < size(); ip++)
    {
        sum_ += sigma_[ip]*operator[](ip);
    }

    return sum_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedFraction<Type, PatchField, GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::fraction() const
{
    const wordList& partNames = PartitionedFraction<Type, PatchField, GeoMesh>::partNames;
    
    // Name of the fraction
    string newName = operator[](0).name();
    const word& partName = partNames[0];
    newName.erase(0, partName.size());
    newName = sigma_[0].name() + newName;
    newName.erase(0, partName.size()+1);
    
    // Geometric field to start from
    GeometricField<Type, PatchField, GeoMesh> fracField
    (
        newName,
        operator[](0)*sigma_[0]
    );
    
    PartitionedFraction<Type, PatchField, GeoMesh> frac
    (
        partNames,
        fracField
    );
    
    for(label ip = 1; ip < size(); ip++)
    {
        frac[ip] = operator[](ip)*sigma_[ip];
    }

    return frac;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).write();
    }

    updateSum();
    PartitionedFraction<Type, PatchField, GeoMesh>::sum_.write();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).readUpdate();
    }
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    Info << "hello from PartitionedField::operator=" << endl;
    PartitionedFraction<Type, PatchField, GeoMesh>(*this) = gf;
    Info << "done" << endl;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator+=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    // Check that you are adding two partitions with the same sigma
    if (&sigma_ != &gf.sigma())
    {
        FatalErrorIn("PartitionedField::operator+=")
            << " attempting to add two fields with different sigmas"
            << abort(FatalError);
    }
    
    PartitionedFraction<Type, PatchField, GeoMesh>::sum_
        += gf.sum();

    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) += gf[ip];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::operator-=
(
    const PartitionedField<Type, PatchField, GeoMesh>& gf
)
{
    // Check that you are subtracting two partitions with the same sigma
    if (&sigma_ != &gf.sigma())
    {
        FatalErrorIn("PartitionedField::operator-=")
            << " attempting to add two fields with different sigmas"
            << abort(FatalError);
    }
    
    PartitionedFraction<Type, PatchField, GeoMesh>::sum_
        -= gf.sum();

    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) -= gf[ip];
    }
}


// ************************************************************************* //
