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
    Foam::limitedCubicUpwindFit

Description
    Upwind biased fit surface interpolation scheme that applies an explicit
    correction to linear.

\*---------------------------------------------------------------------------*/

#ifndef limitedCubicUpwindFit_H
#define limitedCubicUpwindFit_H

#include "UpwindFitData.H"
#include "UpwindFitScheme.H"
#include "linear.H"
#include "cubicUpwindFitPolynomial.H"
#include "upwindCFCCellToFaceStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class limitedCubicUpwindFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Polynomial, class Stencil>
class limitedCubicUpwindFit
:
    public linear<Type>
{
    // Private Data

        //- Reference to the surface flux used to choose upwind direction
        const surfaceScalarField& faceFlux_;

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        limitedCubicUpwindFit(const limitedCubicUpwindFit&);

        //- Disallow default bitwise assignment
        void operator=(const limitedCubicUpwindFit&);


public:

    //- Runtime type information
    TypeName("limitedCubicUpwindFit");


    // Constructors

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        limitedCubicUpwindFit(const fvMesh& mesh, Istream& is)
        :
            linear<Type>(mesh),
            faceFlux_(mesh.lookupObject<surfaceScalarField>(word(is))),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}


        //- Construct from mesh, faceFlux and Istream
        limitedCubicUpwindFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            linear<Type>(mesh),
            faceFlux_(faceFlux),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
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

            const extendedUpwindCellToFaceStencil& stencil = Stencil::New
            (
                mesh,
                false,          //pureUpwind
                scalar(0.5)
            );

            const UpwindFitData<Polynomial>& ufd =
            UpwindFitData<Polynomial>::New
            (
                mesh,
                stencil,
                true,           //calculate as offset to linear
                linearLimitFactor_,
                centralWeight_
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

#define makeUpwindFitSurfaceInterpolationTypeScheme\
(                                                                             \
    SS,                                                                       \
    POLYNOMIAL,                                                               \
    STENCIL,                                                                  \
    TYPE                                                                      \
)                                                                             \
                                                                              \
typedef limitedCubicUpwindFit<TYPE, POLYNOMIAL, STENCIL>                            \
    limitedCubicUpwindFit##TYPE##POLYNOMIAL##STENCIL##_;                            \
defineTemplateTypeNameAndDebugWithName                                        \
    (limitedCubicUpwindFit##TYPE##POLYNOMIAL##STENCIL##_, #SS, 0);                  \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<limitedCubicUpwindFit<TYPE, POLYNOMIAL, STENCIL> >                                 \
    add##SS##STENCIL##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<limitedCubicUpwindFit<TYPE, POLYNOMIAL, STENCIL> >                                 \
    add##SS##STENCIL##TYPE##MeshFluxConstructorToTable_;

#define makeUpwindFitSurfaceInterpolationScheme(SS, POLYNOMIAL, STENCIL)      \
                                                                              \
makeUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,scalar)     \
makeUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,vector)     \
makeUpwindFitSurfaceInterpolationTypeScheme                                   \
(                                                                             \
    SS,                                                                       \
    POLYNOMIAL,                                                               \
    STENCIL,                                                                  \
    sphericalTensor                                                           \
)                                                                             \
makeUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,symmTensor) \
makeUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
