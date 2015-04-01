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
    Foam::meshToPointField

Description
    Templated class for holding the stencils of cells and
    interpolation weights for interpolating from the cell centres of
    a mesh onto a field of points.

SourceFiles
    meshToPointField.C

\*---------------------------------------------------------------------------*/

#ifndef meshToPointField_H
#define meshToPointField_H

#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToPointField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class ApproxType, class ExtendedStencil>
class meshToPointField
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- The stencil the fit is based on
        const ExtendedStencil& cellStencil_;
        
        //- The stencil for each point
        labelListList pointStencil_;

        //- The weights for each cell of the point stencil
        scalarListList weights_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshToPointField(const meshToPointField<Type, ApproxType, ExtendedStencil>&);

        //- Disallow default bitwise assignment
        void operator=(const meshToPointField<Type, ApproxType, ExtendedStencil>&);


public:

    // Constructors

        //- Construct from components
        meshToPointField
        (
            const fvMesh& mesh,
            const pointField& pts,
            const ExtendedStencil& cellStencil__
        );

        //- Construct from mesh (need to set points before using)
        meshToPointField
        (
            const fvMesh& mesh,
            const ExtendedStencil& cellStencil__
        );


    //- Destructor
    ~meshToPointField();


    // Access functions
    
        //- Return reference to the cell stencil
        const ExtendedStencil& cellStencil() const { return cellStencil_; }

        //- Return reference to the point stencil
        const ExtendedStencil& pointStencil() const { return pointStencil_; }

        //- interpolation weights
        const scalarListList& weights() const { return weights_; }

    // Member Functions

        //- Interpolate volField onto these points
        tmp<Field<Type> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
        //- Calculate the stencils and weights for a set of points
        void setPoints(const pointField& pts);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshToPointFieldTemplates.C"
#   include "meshToPointField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
