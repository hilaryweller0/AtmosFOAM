/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::limitedSurfaceInterpolationSchemeDebug

Description
    Abstract base class for limited surface interpolation schemes.

SourceFiles
    limitedSurfaceInterpolationSchemeDebug.C

\*---------------------------------------------------------------------------*/

#ifndef limitedSurfaceInterpolationSchemeDebug_H
#define limitedSurfaceInterpolationSchemeDebug_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class limitedSurfaceInterpolationSchemeDebug Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedSurfaceInterpolationSchemeDebug
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow copy construct
        limitedSurfaceInterpolationSchemeDebug
        (
            const limitedSurfaceInterpolationSchemeDebug&
        );

        //- Disallow default bitwise assignment
        void operator=(const limitedSurfaceInterpolationSchemeDebug&);


protected:

    // Protected data

        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("limitedSurfaceInterpolationSchemeDebug");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationSchemeDebug,
            Mesh,
            (
                const fvMesh& mesh,
                Istream& schemeData
            ),
            (mesh, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationSchemeDebug,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh and faceFlux
        limitedSurfaceInterpolationSchemeDebug
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        limitedSurfaceInterpolationSchemeDebug
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationSchemeDebug<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationSchemeDebug<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~limitedSurfaceInterpolationSchemeDebug();


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        //- Return the interpolation weighting factors for the given field,
        //  by limiting the given weights with the given limiter
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const surfaceScalarField& CDweights,
            tmp<surfaceScalarField> tLimiter
        ) const;

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        flux
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makelimitedSurfaceInterpolationTypeSchemeDebug(SS, Type)                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;                                \
                                                                               \
limitedSurfaceInterpolationSchemeDebug<Type>::addMeshConstructorToTable<SS<Type>>   \
    add##SS##Type##MeshConstructorToLimitedTable_;                             \
                                                                               \
limitedSurfaceInterpolationSchemeDebug<Type>::                                      \
    addMeshFluxConstructorToTable<SS<Type>>                                    \
    add##SS##Type##MeshFluxConstructorToLimitedTable_;

#define makelimitedSurfaceInterpolationSchemeDebug(SS)                              \
                                                                               \
makelimitedSurfaceInterpolationTypeSchemeDebug(SS, scalar)                          \
makelimitedSurfaceInterpolationTypeSchemeDebug(SS, vector)                          \
makelimitedSurfaceInterpolationTypeSchemeDebug(SS, sphericalTensor)                 \
makelimitedSurfaceInterpolationTypeSchemeDebug(SS, symmTensor)                      \
makelimitedSurfaceInterpolationTypeSchemeDebug(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "limitedSurfaceInterpolationSchemeDebug.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
