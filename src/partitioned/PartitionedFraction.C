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
#include "PartitionedField.H"

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
            "sum." + io.name(),
            io.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            io.writeOpt()
        ),
        mesh
    )
{
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    partNames[ip] + '.' + io.name(),
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
    sum_.rename("sum."+field.name());

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
                    partNames[ip] + '.' + field.name(),
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
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::PartitionedFraction
(
    const PartitionedFraction<Type, PatchField, GeoMesh>& copyFrom
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(copyFrom.size()),
    partNames(copyFrom.partNames),
    sum_(copyFrom.sum())
{
    // Set fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                copyFrom[ip]
            )
        );
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::PartitionedFraction
(
    const word& newName,
    const PartitionedFraction<Type, PatchField, GeoMesh>& copyFrom
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(copyFrom.size()),
    partNames(copyFrom.partNames),
    sum_(copyFrom.sum())
{
    // Set fields for each partition
    for(label ip = 0; ip < size(); ip++)
    {
        (*this).set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                partNames[ip] + '.' + newName,
                copyFrom[ip]
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
Foam::PartitionedField<Type, PatchField, GeoMesh>
Foam::PartitionedFraction<Type, PatchField, GeoMesh>::field
(
    const PartitionedFraction<Type, PatchField, GeoMesh>& sigma
) const
{
    // Name of the new field
    string newName = operator[](0).name();
    newName.erase(0, sigma[0].name().size()+1);

    // Geometric field to start from
    GeometricField<Type, PatchField, GeoMesh> ff
    (
        newName,
        operator[](0)/sigma[0]
    );
    
    PartitionedField<Type, PatchField, GeoMesh> f(partNames, ff, sigma);
    
    for(label ip = 1; ip < size(); ip++)
    {
        Info << "Setting partition " << ip << endl;
        f[ip] = operator[](ip)/sigma[ip];
    }
    f.updateSum();

    return f;
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


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedFraction<Type, PatchField, GeoMesh>::operator=
(
    const PartitionedFraction<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    Info << "hello from PartitionedFraction::operator=" << endl;
    PtrList<GeometricField<Type, PatchField, GeoMesh>>(*this) = gf;
    Info << "About to update sum" << endl;
    updateSum();
    Info << "Done" << endl;
}

// ************************************************************************* //
