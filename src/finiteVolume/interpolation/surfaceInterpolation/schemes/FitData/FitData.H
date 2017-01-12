/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::FitData

Description
    Data for the upwinded and centred polynomial fit interpolation schemes.
    The linearCorrection_ determines whether the fit is for a corrected
    linear scheme (first two coefficients are corrections for owner and
    neighbour) or a pure upwind scheme (first coefficient is correction for
    owner; weight on face taken as 1).

SourceFiles
    FitData.C

\*---------------------------------------------------------------------------*/

#ifndef FitData_H
#define FitData_H

#include "autoPtr.H"
#include "MeshObject.H"
#include "fvMesh.H"
#include "fitResult.H"
#include "fitWeights.H"
#include "fitCoefficients.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FitData Declaration
\*---------------------------------------------------------------------------*/

template<class FitDataType, class ExtendedStencil, class Polynomial>
class FitData
:
    public MeshObject<fvMesh, MoveableMeshObject, FitDataType>
{
    // Private data

        //- The stencil the fit is based on
        const ExtendedStencil& stencil_;

        //- Is scheme correction on linear (true) or on upwind (false)
        const bool linearCorrection_;

        //- Weights for central stencil
        const scalar centralWeight_;

        //- Dimensionality of the geometry
        const label dim_;

        //- Determines the orientation of the local stencil coordinates
        //- For spherical geometry, a plane tangent to the face centre is used
        //- Otherwise, the jdir is chosen to be maximimally orthogonal to idir
        const bool sphericalGeometry_;


protected:

        //- Find the normal direction (i) and j and k directions for face faci
        void findFaceDirs
        (
            vector& idir,        // value changed in return
            vector& jdir,        // value changed in return
            vector& kdir,        // value changed in return
            const label faci,
            const List<point>& C
        );

public:

    // Constructors

        //- Construct from components
        FitData
        (
            const fvMesh& mesh,
            const ExtendedStencil& stencil,
            const bool linearCorrection,
            const scalar centralWeight,
            const bool sphericalGeometry = false
        );


    //- Destructor
    virtual ~FitData()
    {}


    // Member functions

        //- Return reference to the stencil
        const ExtendedStencil& stencil() const
        {
            return stencil_;
        }

        //- Return weight for central stencil
        scalar centralWeight() const
        {
            return centralWeight_;
        }

        //- Dimensionality of the geometry
        label dim() const
        {
            return dim_;
        }

        bool linearCorrection() const
        {
            return linearCorrection_;
        }

        autoPtr<fitResult> calcFit
        (
            fitCoefficients& coefficients, // coefficients to be set
            const List<point>& stencilPoints,  // Stencil points
            const label faci,     // Current face index
            bool owner
        );

        //- Calculate the fit for the specified face and set the coefficients
        void calcFit
        (
            scalarList& coeffsi, // coefficients to be set
            scalarList& wts,
            const List<point>&,  // Stencil points
            const scalar wLin,   // Weight for linear approximation (weights
                                 // nearest neighbours)
            const label faci    // Current face index
        );

        //- Calculate the fit for all the faces
        virtual void calcFit() = 0;

        //- Recalculate weights (but not stencil) when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
