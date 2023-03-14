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
    Foam::faceToPoint

Description
    Abstract base type for holding the interpolation weights for interpolating
    from the face centres of a mesh onto the points of the same mesh

SourceFiles
    faceToPoint.C
    faceToPointI.H

\*---------------------------------------------------------------------------*/

#ifndef faceToPoint_H
#define faceToPoint_H

#include "meshToPoint.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faceToPoint Declaration
\*---------------------------------------------------------------------------*/

template<class approxType>
class faceToPoint
{
    // Private data

        const fvMesh& mesh_;
        
        scalarListList weights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceToPoint(const faceToPoint<approxType>&);

        //- Disallow default bitwise assignment
        void operator=(const faceToPoint<approxType>&);


public:

    // Constructors

        //- Construct from components
        faceToPoint(const fvMesh& mesh);


    //- Destructor
    ~faceToPoint();


    // Access functions
    
        //- interpolation weights
        const scalarListList& weights() const { return weights_; }

    // Member Functions

        //- Interpolate surfaceField onto these points
        template<class Type>
        inline Field<Type> interpolate
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceToPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faceToPoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
