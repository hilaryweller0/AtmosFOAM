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
    A list of GeometricFields, one for each partion of the atmosphere

SourceFiles
    PartitionedField.C

\*---------------------------------------------------------------------------*/

#ifndef PartitionedField_H
#define PartitionedField_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PartitionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class PartitionedField
:
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Public data
public:
        //- The names of the parts
        const wordList& partNames;


    // Private data
private:
        //- The volume fractions in each part
        PartitionedField<Type, PatchField, GeoMesh>& sigma_;

        //- The sum over all partitions
        GeometricField<Type, PatchField, GeoMesh> sum_;
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        PartitionedField(const PartitionedField&);

        //- Disallow default bitwise assignment
        void operator=(const PartitionedField&);


public:

    // Public typedefs

        //- Type of mesh on which this PartitionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

    // Constructors

        //- Construct from partion names and IOobject (must read in fields)
        PartitionedField
        (
            const wordList& partNames__,
            const IOobject& io,
            const Mesh& mesh
        );

        //- Construct from partion names and the volume fraction (sigma)
        //- (must read in the fields)
        PartitionedField
        (
            const wordList& partNames__,
            const IOobject& io,
            const Mesh& mesh,
            const PartitionedField& sigma__
        );

        //- Construct from partion names and IOobject.
        //- Initialise all parts to "field"
        PartitionedField
        (
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field
        );

        //- Construct from partion names and the volume fraction (sigma)
        //- Initialise all parts to "field"
        PartitionedField
        (
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field,
            const PartitionedField& sigma__
        );


    //- Destructor
    ~PartitionedField();

    // Access functions
    
        //- Sigma - volume fractions in each partition
        const PartitionedField<Type, PatchField, GeoMesh>& sigma() const
            { return sigma_; }

        //- Return the sum without updating
        const GeometricField<Type, PatchField, GeoMesh>& sum() const
            { return sum_; }

    // Member Functions

        //- Sum over partitions and return the sum
        const GeometricField<Type, PatchField, GeoMesh>& updateSum();

        // Output
        void write();

        // Input (re-read all fields)
        void readUpdate();
        
        // Functions from PtrList
        label size() { return (*this).size();}
        
        const GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label ip
        )
        { return (*this).operator[](ip);}
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
