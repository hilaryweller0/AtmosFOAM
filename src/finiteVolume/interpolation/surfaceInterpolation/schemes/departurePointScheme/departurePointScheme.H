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
    Foam::departurePointScheme

Description
    Interpolates onto a point upstream of the face - an explicit correction
    on linear

\*---------------------------------------------------------------------------*/

#ifndef departurePointScheme_H
#define departurePointScheme_H

#include "departurePointData.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class departurePointScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class ApproxType, class Stencil>
class departurePointScheme
:
    public linear<Type>
{
    // Private Data

        //- Reference to the departure points
        const surfaceVectorField& departurePoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        departurePointScheme(const departurePointScheme&);

        //- Disallow default bitwise assignment
        void operator=(const departurePointScheme&);


public:

    //- Runtime type information
    TypeName("departurePointScheme");


    // Constructors

        //- Construct from mesh and Istream
        //  The name of the departure points is read from the 
        //  Istream and looked-up from the mesh objectRegistry
        departurePointScheme(const fvMesh& mesh, Istream& is)
        :
            linear<Type>(mesh),
            departurePoints_(mesh.lookupObject<surfaceVectorField>(word(is)))
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            const extendedCentredCellToCellStencil& stencil = Stencil::New
            (
                mesh
            );

            const departurePointData<ApproxType>& ufd
                = departurePointData<ApproxType>::New
            (
                mesh,
                stencil,
                departurePoints_
            );

            return ufd.interpolate(vf) - linearInterpolate(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeDeparturePointSurfaceInterpolationTypeScheme\
(                                                                             \
    SS,                                                                       \
    APPROXTYPE,                                                               \
    STENCIL,                                                                  \
    TYPE                                                                      \
)                                                                             \
                                                                              \
typedef departurePointScheme<TYPE, APPROXTYPE, STENCIL>                            \
    departurePointScheme##TYPE##APPROXTYPE##STENCIL##_;                            \
defineTemplateTypeNameAndDebugWithName                                        \
    (departurePointScheme##TYPE##APPROXTYPE##STENCIL##_, #SS, 0);                  \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<departurePointScheme<TYPE, APPROXTYPE, STENCIL> >                                 \
    add##SS##STENCIL##TYPE##MeshConstructorToTable_;

#define makeDeparturePointSurfaceInterpolationScheme(SS, APPROXTYPE, STENCIL)      \
                                                                              \
makeDeparturePointSurfaceInterpolationTypeScheme(SS,APPROXTYPE,STENCIL,scalar)     \
makeDeparturePointSurfaceInterpolationTypeScheme(SS,APPROXTYPE,STENCIL,vector)     \
makeDeparturePointSurfaceInterpolationTypeScheme                                   \
(                                                                             \
    SS,                                                                       \
    APPROXTYPE,                                                               \
    STENCIL,                                                                  \
    sphericalTensor                                                           \
)                                                                             \
makeDeparturePointSurfaceInterpolationTypeScheme(SS,APPROXTYPE,STENCIL,symmTensor) \
makeDeparturePointSurfaceInterpolationTypeScheme(SS,APPROXTYPE,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
