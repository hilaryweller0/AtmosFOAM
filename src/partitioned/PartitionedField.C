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
    const word& baseName__,
    const wordList& partNames__,
    const Mesh& mesh,
    const word& timeName
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            "sum."+baseName_, timeName, mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensioned<Type>("sum", dimless, pTraits<Type>::zero)
    ),
    needsSigma_(false),
    sigma_(*this),
    ddtPtr_(NULL)
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
                    partNames()[ip] + '.' + baseName(), timeName, mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const Mesh& mesh,
    const word& timeName,
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            "sum."+baseName_, timeName, mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensioned<Type>("sum", dimless, pTraits<Type>::zero)
    ),
    needsSigma_(true),
    sigma_(sigma__),
    ddtPtr_(NULL)
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
                    partNames()[ip] + '.' + baseName(), timeName, mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    updateSum();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            "sum."+baseName_, field.time().timeName(), field.mesh(),
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        field.mesh(),
        dimensioned<Type>("sum", dimless, pTraits<Type>::zero)
    ),
    needsSigma_(false),
    sigma_(*this),
    ddtPtr_(NULL)
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
                    partNames()[ip] + '.' + baseName(), 
                    field.mesh().time().timeName(),
                    field.mesh(),
                    IOobject::NO_READ, IOobject::AUTO_WRITE
                ),
                field
            )
        );
    }
    updateSum();
}



template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::PartitionedField
(
    const word& baseName__,
    const wordList& partNames__,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__
)
:
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(partNames__.size()),
    baseName_(baseName__),
    partNames_(partNames__),
    sum_
    (
        IOobject
        (
            "sum."+baseName_, field.time().timeName(), field.mesh(),
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        field.mesh(),
        dimensioned<Type>("sum", dimless, pTraits<Type>::zero)
    ),
    needsSigma_(true),
    sigma_(sigma__),
    ddtPtr_(NULL)
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
                    partNames_[ip]+'.'+baseName_,
                    field.mesh().time().timeName(),
                    field.mesh(),
                    IOobject::NO_READ, IOobject::AUTO_WRITE
                ),
                field
                
            )
        );
    }
    updateSum();
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
    // Sum depends on whether it is sigma weighted
    
    if (!needsSigma_)
    {
        sum_.dimensions().reset(operator[](0).dimensions());
        sum_ = operator[](0);

        // Sum contributions from other partitions
        for (label ip = 1; ip < size(); ip++)
        {
            sum_ += operator[](ip);
        }
    }
    else
    {
        sum_.dimensions().reset
        (
            operator[](0).dimensions()*sigma_[0].dimensions()
        );
    
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
Foam::PartitionedField<Type, PatchField, GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::timesSigma() const
{
    if (!needsSigma())
    {
        FatalErrorIn("PartitionedField::timesSigma")
            << "cannot calculate timesSigma without sigma"
            << abort(FatalError);
    }
    
    // Name of the timesSigma
    string newName = sigma_.baseName()+'.'+baseName();
    
    // Geometric field to start from
    const GeometricField<Type, PatchField, GeoMesh>& f = operator[](0);
    GeometricField<Type, PatchField, GeoMesh> fracField
    (
        IOobject(newName, f.time().timeName(), f.mesh()),
        f*sigma_[0],
        f.boundaryField().types()
    );
    
    PartitionedField<Type, PatchField, GeoMesh> frac
    (
        newName, partNames(), fracField
    );
    
    for(label ip = 1; ip < size(); ip++)
    {
        frac[ip] = operator[](ip)*sigma_[ip];
    }
    
    frac.updateSum();
    
    return frac;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>
Foam::PartitionedField<Type, PatchField, GeoMesh>::divideBy
(
    const PartitionedField<scalar, PatchField, GeoMesh>& sigma__,
    const word newBaseName
) const
{
    if (needsSigma())
    {
        FatalErrorIn("PartitionedField::divideBy")
            << "cannot calculate divideBy sigma for a field with sigma"
            << abort(FatalError);
    }

    // Name of the new field
    string newName = newBaseName;
    if (newBaseName == "")
    {
        newName = baseName();
        newName.erase(0, sigma()[0].name().size()+1);
    }
    
    // Geometric field to start from
    GeometricField<Type, PatchField, GeoMesh> ff
    (
        newName, operator[](0)/sigma__[0]
    );

    // New partitionedField
    PartitionedField<Type, PatchField, GeoMesh> f
    (
        newName, partNames(), ff, sigma__
    );
    
    for(label ip = 1; ip < size(); ip++)
    {
        f[ip] = operator[](ip)/sigma__[ip];
    }

    return f;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::PartitionedField<Type, PatchField, GeoMesh>::storeTime()
{
    // Only call this function once (only set the pointers once)
    if (ddtPtr_)
    {
        FatalErrorIn("PartitionedField::storeTime")
            << "already called, ddtPtr_ already set"
            << abort(FatalError);
    }

    // Store old time for all of the fields
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip).oldTime();
    }
    
    // Store old sum
    sum().oldTime();
    
    // Create fields for all the time directories
    ddtPtr_ = new PtrList<GeometricField<Type, PatchField, GeoMesh> >
    (
        partNames().size()
    );
    for(label ip = 0; ip < size(); ip++)
    {
        ddt().set
        (
            ip,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    partNames()[ip]+".ddt."+baseName(),
                    operator[](ip).mesh().time().timeName(),
                    operator[](ip).mesh()
                ),
                (operator[](ip) - operator[](ip).oldTime())
                    /operator[](ip).time().deltaT(),
                "fixedValue"
            )
        );
        ddt()[ip].oldTime();
    }
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

    // Write out old time and rate of change if set
    if (ddtPtr_)
    {
        for(label ip = 0; ip < size(); ip++)
        {
            operator[](ip).oldTime().write();
            ddt()[ip].write();
        }
    }
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
    
    //PtrList<GeometricField<Type, PatchField, GeoMesh>>(*this) = gf;
    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) = gf[ip];
    }
    sum_ = gf.sum();
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
    
    sum_ += gf.sum();

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
    
    sum_ -= gf.sum();

    for(label ip = 0; ip < size(); ip++)
    {
        operator[](ip) -= gf[ip];
    }
}


// ************************************************************************* //
