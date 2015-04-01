/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::approxType

Description
    abstract base type for approximation types such as polynomial fits and
    radial basis functions

SourceFiles

    approxType.C

\*---------------------------------------------------------------------------*/

#ifndef approxType_H
#define approxType_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class approxType Declaration
\*---------------------------------------------------------------------------*/

class approxType
{
public:

    //- Destructor
    virtual ~approxType(){}

    // Member functions

        //- Find local coordinate directions for a given stencil onto a given point
        void findDirs
        (
            List<vector>& dir, const pointField& points,
            const labelList& stencil, const point& target, const int nDims
        ) const;

        //- Calculate the approximattion weights from cell centres.
        //- If centralCell is not -1, this cell is weighted in the fit, 
        //- along with the zeroth cell of the stencil
        virtual int calcWeights
        (
            scalarList& weights,
            const polyMesh& mesh,
            const labelList& stencil,
            const point& target,
            const int centralCell = -1
        ) const = 0;

//        //- Calculate the approximattion weights using given direction
//        virtual int calcWeights
//        (
//            scalarList& weights, const polyMesh& mesh,
//            const labelList& stencil, const point& p,
//            const List<vector>& dir
//        ) const = 0;

        //- Calculate the approximattion weights from pointField
        virtual int calcWeights
        (
            scalarList& weights, const pointField& points,
            const labelList& stencil, const point& target,
            const int nDims, const int centralCell = -1
        ) const = 0;

        //- Calculate the approximattion weights from pointField using all points
        virtual int calcWeights
        (
            scalarList& weights, const pointField& points,
            const point& target, const int nDims, const int centralCell = -1
        ) const = 0;

//        //- Calculate the weights for interpolating from the mesh cell centres
//        //- using given direction
//        virtual int calcWeights
//        (
//            scalarList& weights, const pointField& points,
//            const labelList& stencil, const point& target,
//            const List<vector>& dir, const int nDims
//        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
