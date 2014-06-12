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
    Foam::meshToPoint

Description
    Abstract base type for holding the stencil of cells and
    interpolation weights for interpolating from the cell centres of
    a mesh onto the given point

SourceFiles
    meshToPoint.C
    meshToPoinI.H

\*---------------------------------------------------------------------------*/

#ifndef meshToPoint_H
#define meshToPoint_H

#include "fvMesh.H"
#include "volFields.H"
#include "approxType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToPoint Declaration
\*---------------------------------------------------------------------------*/

template<class approxType>
class meshToPoint
:
    public point
{
    // Private data

        const fvMesh& mesh_;
        
        labelList stencil_;
        
        scalarList weights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshToPoint(const meshToPoint<approxType>&);

        //- Disallow default bitwise assignment
        void operator=(const meshToPoint<approxType>&);


public:

    // Constructors

        //- Construct from components
        meshToPoint(const point& p, const fvMesh& mesh);


    //- Destructor
    ~meshToPoint();


    // Access functions
    
        //- interpolation stencils
        const labelList& stencil() const { return stencil_; }

        //- interpolation weights
        const scalarList& weights() const { return weights_; }

    // Member Functions

        //- Interpolate volField onto this point
        template<class Type>
        inline Type interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshToPoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
