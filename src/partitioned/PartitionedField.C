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
//#include "moreListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const IOobject& io,
    const Mesh& mesh
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames.size()),
    partNames(partNames__),
    sigma_(*this),
    // Read in the sum
    sum_
    (
        IOobject
        (
            io.name()+"Sum",
            io.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            io.writeOpt()
        ),
        mesh
    )
{
    // Read in fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    io.name()+partNames[ip],
                    io.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    io.writeOpt()
                ),
                mesh
            )
        );
    }
    
    // Update the sum
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const IOobject& io,
    const Mesh& mesh,
    const PartitionedField& sigma__
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames.size()),
    partNames(partNames__),
    sigma_(sigma__),
    // Read in the sum
    sum_
    (
        IOobject
        (
            io.name()+"Sum",
            io.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            io.writeOpt()
        ),
        mesh
    )
{
    // Read in fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    io.name()+partNames[ip],
                    io.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    io.writeOpt()
                ),
                mesh
            )
        );
    }
    
    // Update the sum
    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames.size()),
    partNames(partNames__),
    sigma_(*this),
    sum_(field)
{
    // Set fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    field.name()+partNames[ip],
                    field.time().timeName(),
                    field.mesh(),
                    IOobject::NO_READ
                ),
                field
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const PartitionedField& sigma__
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames.size()),
    partNames(partNames__),
    sigma_(sigma__),
    sum_(field)
{
    // Set fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    field.name()+partNames[ip],
                    field.time().timeName(),
                    field.mesh(),
                    IOobject::NO_READ
                ),
                field
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::~PartitionedField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::PartitionedField<Type, PatchField, GeoMesh>::updateSum()
{
    // The sum depends on whether the partionedField is sigma
    if (this == &sigma_)
    {
        sum_ = sigma_[0];

        // Sum contributions from other partitions
        for (label ip = 1; ip < size(); ip++)
        {
            sum_ += sigma_[ip];
        }
    }
    else
    {
        sum_ = sigma_[0]*operator[](0);

        // Sum contributions from other partitions
        for (label ip = 1; ip < size(); ip++)
        {
            sum_ += sigma_[ip]*operator[](ip);
        }
    }
    return sum_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).write();
    }

    updateSum();
    sum_.write();
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


// ************************************************************************* //
