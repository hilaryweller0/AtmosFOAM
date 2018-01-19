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

Class
    Foam::PartitionedField

Description
    A list of GeometricFields, one for each partion of the atmosphere. A
    PartitionedField is a PartitioneFraction plus a sigma field which gives
    the fraction so that
    sum = sum(sigma*field)

SourceFiles
    PartitionedField.C

\*---------------------------------------------------------------------------*/

#ifndef PartitionedField_H
#define PartitionedField_H

#include "PartitionedFraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PartitionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class PartitionedField
:
    public PartitionedFraction<Type, PatchField, GeoMesh>
{
    // Private data

        //- The fractions in each partition
        const PartitionedFraction<Type, PatchField, GeoMesh>& sigma_;


public:

    // Public typedefs

        //- Type of mesh on which this PartitionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;
        
    // Constructors

        //- Construct from partion names and IOobject (must read in fields)
        //- and sigma
        PartitionedField
        (
            const wordList& partNames__,
            const IOobject& io,
            const Mesh& mesh,
            const PartitionedFraction<Type, PatchField, GeoMesh>& sigma__
        );

        //- Construct from partion names and IOobject and sigma.
        //- Initialise all parts to "field"
        PartitionedField
        (
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field,
            const PartitionedFraction<Type, PatchField, GeoMesh>& sigma__
        );

        //- Copy construct
        PartitionedField(const PartitionedField&);

        //- Copy and rename
        PartitionedField
        (
            const word& newName,
            const PartitionedField<Type, PatchField, GeoMesh>& copyFrom
        );

    //- Destructor
    ~PartitionedField();

    // Access functions
    
        //- Return sigma
        const PartitionedFraction<Type, PatchField, GeoMesh>& sigma() const
            { return sigma_; }

    // Member Functions

        //- Sum over partitions weighting by sigma and return the sum
        const GeometricField<Type, PatchField, GeoMesh>& updateSum();

        //- Multiply the field by sigma
        PartitionedFraction<Type, PatchField, GeoMesh> fraction() const;

        // Output
        void write();

        // Input (re-read all fields)
        void readUpdate();

    // Member operators
    
        void operator=(const PartitionedField<Type, PatchField, GeoMesh>&);
        void operator+=(const PartitionedField<Type, PatchField, GeoMesh>&);
        void operator-=(const PartitionedField<Type, PatchField, GeoMesh>&);

    // Functions from PtrList

        label size() const
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::size();
        }
        
        const GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label ip
        ) const
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](ip);
        }

        GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label ip
        )
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](ip);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PartitionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
