/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::departurePointData

Description
    Data for the departure point scheme for interpolating from volFields
    onto the departure points of faces

SourceFiles
    departurePointData.C

\*---------------------------------------------------------------------------*/

#ifndef departurePointData_H
#define departurePointData_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedCentredCellToCellStencil;

/*---------------------------------------------------------------------------*\
                    Class departurePointData Declaration
\*---------------------------------------------------------------------------*/

template<class ApproxType>
class departurePointData
:
    public MeshObject<fvMesh, MoveableMeshObject, departurePointData<ApproxType> >
{
    // Private data

        //- The stencil the fit is based on
        const extendedCentredCellToCellStencil& stencil_;

        //- Reference to the departure points
        const surfaceVectorField& departurePoints_;

        //- The cell that each departure point is in (to find the stencil)
        labelList pointInCell_;

        //- For each departure point, the coefficients to multiply the
        //  stencil cell values by
        List<scalarList> coeffs_;


public:

    TypeName("departurePointData");


    // Constructors

        //- Construct from components
        departurePointData
        (
            const fvMesh& mesh,
            const extendedCentredCellToCellStencil& stencil,
            const surfaceVectorField& depPoints
        );


    //- Destructor
    virtual ~departurePointData()
    {}


    // Member functions

        //- Calculate the stencils and weights for a set of points
        void setPoints();

        //- Interpolate volField onto these points
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
        //- Recalculate weights (but not stencil) when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "departurePointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
