/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::meshToPoints

Description
    Abstract base type for holding the interpolation weights for interpolating
    from the cell centres of a mesh onto the points of the same mesh

SourceFiles
    meshToPoints.C
    meshToPoinI.H

\*---------------------------------------------------------------------------*/

#ifndef meshToPoints_H
#define meshToPoints_H

#include "meshToPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToPoints Declaration
\*---------------------------------------------------------------------------*/

template<class approxType>
class meshToPoints
{
    // Private data

        const fvMesh& mesh_;
        
        scalarListList weights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshToPoints(const meshToPoints<approxType>&);

        //- Disallow default bitwise assignment
        void operator=(const meshToPoints<approxType>&);


public:

    // Constructors

        //- Construct from components
        meshToPoints(const fvMesh& mesh);


    //- Destructor
    ~meshToPoints();


    // Access functions
    
        //- interpolation weights
        const scalarListList& weights() const { return weights_; }

    // Member Functions

        //- Interpolate volField onto these points
        template<class Type>
        inline Field<Type> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshToPoints.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
