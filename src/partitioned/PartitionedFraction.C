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

#include "PartitionedFraction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::PartitionedFraction
(
    const wordList& partNames__,
    const IOobject& io,
    const Mesh& mesh
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    partNames(partNames__),
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
    Info << "Reading in " << io.name() << " for each partition" << endl;
    for(label ip = 0; ip < size(); ip++)
    {
        Info << "Partition " << ip << " called " << partNames[ip] << endl;
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

    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::PartitionedFraction
(
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    partNames(partNames__),
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
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::~PartitionedFraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::updateSum()
{
    sum_ = operator[](0);

    // Sum contributions from other partitions
    for (label ip = 1; ip < size(); ip++)
    {
        sum_ += operator[](ip);
    }

    return sum_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedFraction<Type, PatchField, GeoMesh>::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).write();
    }

    updateSum();
    sum_.write();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedFraction<Type, PatchField, GeoMesh>::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).readUpdate();
    }
    updateSum();
}


// ************************************************************************* //
