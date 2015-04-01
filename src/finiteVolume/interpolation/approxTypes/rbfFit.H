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
    Foam::rbfFit

Description
    Calculates the weights for a radial basis function fit given a stencil of
    points within a mesh and a point to interpolate onto

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef rbfFit_H
#define rbfFit_H

#include "approxType.H"
#include "ListOps.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rbfFit Declaration
\*---------------------------------------------------------------------------*/

class rbfFit
:
    public approxType
{
private:
    // Private member function
        //- The radial basis function
        scalar rbf(const scalar r) const { return sqrt(sqr(r) + 1); }

public:
        //- Destructor
        virtual ~rbfFit(){}


        //- Member functions

        //- Calculate the approximattion weights from cell centres
        int calcWeights
        (
            scalarList& weights,
            const polyMesh& mesh,
            const labelList& stencil,
            const point& target,
            const int centralCell = -1
        ) const;

        //- Calculate the approximattion weights from pointField using all points
        int calcWeights
        (
            scalarList& weights, const pointField& points,
            const point& target, const int nDims, const int centralCell = -1
        ) const;

        //- Calculate the approximattion weights from pointField
        int calcWeights
        (
            scalarList& weights, const pointField& points,
            const labelList& stencil, const point& target,
            const int nDims, const int centralCell = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
