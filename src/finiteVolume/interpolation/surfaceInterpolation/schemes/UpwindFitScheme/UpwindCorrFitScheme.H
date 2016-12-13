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
    Foam::UpwindCorrFitScheme

Description
    Upwind biased fit surface interpolation scheme that applies an explicit
    correction to upwind.

\*---------------------------------------------------------------------------*/

#ifndef UpwindCorrFitScheme_H
#define UpwindCorrFitScheme_H

#include "UpwindCorrFitData.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UpwindCorrFitScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Polynomial, class Stencil>
class UpwindCorrFitScheme
:
    public upwind<Type>
{
    // Private Data

        //- Reference to the surface flux used to choose upwind direction
        const surfaceScalarField& faceFlux_;

        //- Weights for central stencil
        const scalar centralWeight_;

        const bool sphericalGeometry_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        UpwindCorrFitScheme(const UpwindCorrFitScheme&);

        //- Disallow default bitwise assignment
        void operator=(const UpwindCorrFitScheme&);


public:

    //- Runtime type information
    TypeName("UpwindCorrFitScheme");


    // Constructors

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        UpwindCorrFitScheme(const fvMesh& mesh, Istream& is)
        :
            upwind<Type>
            (
                mesh,
                mesh.lookupObject<surfaceScalarField>(word(is))
            ),
            faceFlux_(upwind<Type>::faceFlux_),
            centralWeight_(1000),
            sphericalGeometry_(readBool(is))
        {}


        //- Construct from mesh, faceFlux and Istream
        UpwindCorrFitScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux),
            faceFlux_(faceFlux),
            centralWeight_(1000),
            sphericalGeometry_(readBool(is))
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

            const extendedUpwindCellToFaceStencilNew& stencil = Stencil::New
            (
                mesh,
                false,          //not pureUpwind
                scalar(0.5),
                false           // not a correction on linear
            );

            const UpwindCorrFitData<Polynomial>& ufd =
            UpwindCorrFitData<Polynomial>::New
            (
                mesh,
                stencil,
                false,           // not calculated as offset to linear
                centralWeight_,
                sphericalGeometry_
            );

            const List<scalarList>& fo = ufd.owncoeffs();
            const List<scalarList>& fn = ufd.neicoeffs();

            return stencil.weightedSum(faceFlux_, vf, fo, fn);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeUpwindCorrFitSurfaceInterpolationTypeScheme\
(                                                                             \
    SS,                                                                       \
    POLYNOMIAL,                                                               \
    STENCIL,                                                                  \
    TYPE                                                                      \
)                                                                             \
                                                                              \
typedef UpwindCorrFitScheme<TYPE, POLYNOMIAL, STENCIL>                            \
    UpwindCorrFitScheme##TYPE##POLYNOMIAL##STENCIL##_;                            \
defineTemplateTypeNameAndDebugWithName                                        \
    (UpwindCorrFitScheme##TYPE##POLYNOMIAL##STENCIL##_, #SS, 0);                  \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<UpwindCorrFitScheme<TYPE, POLYNOMIAL, STENCIL> >                                 \
    add##SS##STENCIL##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<UpwindCorrFitScheme<TYPE, POLYNOMIAL, STENCIL> >                                 \
    add##SS##STENCIL##TYPE##MeshFluxConstructorToTable_;

#define makeUpwindCorrFitSurfaceInterpolationScheme(SS, POLYNOMIAL, STENCIL)      \
                                                                              \
makeUpwindCorrFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,scalar)     \
makeUpwindCorrFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,vector)     \
makeUpwindCorrFitSurfaceInterpolationTypeScheme                                   \
(                                                                             \
    SS,                                                                       \
    POLYNOMIAL,                                                               \
    STENCIL,                                                                  \
    sphericalTensor                                                           \
)                                                                             \
makeUpwindCorrFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,symmTensor) \
makeUpwindCorrFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
