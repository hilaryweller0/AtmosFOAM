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
    Foam::polyFit

Description
    Calculates the weights for a polynomial fit given a stencil of points
    withing a mesh and a point to interpolate onto

SourceFiles

    polyFit.C

\*---------------------------------------------------------------------------*/

#ifndef polyFit_H
#define polyFit_H

#include "approxType.H"
#include "ListOps.H"
#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyFit Declaration
\*---------------------------------------------------------------------------*/

//- Enumerations for the order (since you can have intermediate orders)
enum orderType{oZERO, oONE, oONEPLUS, oTWO, oTWOPLUS, oTWOPLUSPLUS, oTHREE, oTHREEPLUS};

template<orderType Order>
class polyFit
:
    public approxType
{
private:
    // Private member function
        //- The minimum size for an dim dimensional polynomial of this order
        //- or of order ord
        int minPolySize(const int dim, const orderType ord = Order) const;


public:

        //- Destructor
        virtual ~polyFit(){}

        //- Member functions

        //- Calculate the approximattion weights from cell centres
        int calcWeights
        (
            scalarList& weights,
            const polyMesh& mesh,
            const labelList& stencil,
            const point& target,
            const int centralCell = 0
        ) const;

        //- Calculate the approximattion weights using given direction
        int calcWeights
        (
            scalarList& weights, const polyMesh& mesh,
            const labelList& stencil, const point& p,
            const List<vector>& dir,
            const int centralCell = 0
        ) const;

        //- Calculate the approximattion weights from pointField
        int calcWeights
        (
            scalarList& weights, const pointField& points,
            const labelList& stencil, const point& target,
            const int nDims,
            const int centralCell = 0
        ) const;

        //- Calculate the approximattion weights from pointField using all points
        int calcWeights
        (
            scalarList& weights, const pointField& points,
            const point& target, const int nDims, const int centralCell = 0
        ) const;

        //- Calculate the weights for interpolating from the mesh cell centres
        //- using given direction
        int calcWeights
        (
            scalarList& weights, const pointField& points,
            const labelList& stencil, const point& target,
            const List<vector>& dir, const int nDims,
            const int centralCell = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "polyFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
