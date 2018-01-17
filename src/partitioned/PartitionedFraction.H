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
    Foam::PartitionedFraction

Description
    A list of GeometricFields, one for each partion of the atmosphere. The
    Fields are fractions of the atmosphere and can be summed to the total

SourceFiles
    PartitionedFraction.C

\*---------------------------------------------------------------------------*/

#ifndef PartitionedFraction_H
#define PartitionedFraction_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PartitionedFraction Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class PartitionedFraction
:
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Public data
public:
        //- The names of the parts
        const wordList& partNames;


    // Protected data
protected:

        //- The sum over all partitions
        GeometricField<Type, PatchField, GeoMesh> sum_;
 
    // Private data
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PartitionedFraction(const PartitionedFraction&);

        //- Disallow default bitwise assignment
        void operator=(const PartitionedFraction&);


public:

    // Public typedefs

        //- Type of mesh on which this PartitionedFraction is instantiated
        typedef typename GeoMesh::Mesh Mesh;

    // Constructors

        //- Construct from partion names and IOobject (must read in fields)
        PartitionedFraction
        (
            const wordList& partNames__,
            const IOobject& io,
            const Mesh& mesh
        );

        //- Construct from partion names and IOobject.
        //- Initialise all parts to "field"
        PartitionedFraction
        (
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field
        );


    //- Destructor
    ~PartitionedFraction();

    // Access functions
    
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PartitionedFraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
