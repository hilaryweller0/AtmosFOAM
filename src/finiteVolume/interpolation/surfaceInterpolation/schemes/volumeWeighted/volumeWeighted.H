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
    Foam::volumeWeighted

Description
    Mid-point interpolation (weighting factors = 0.5) scheme class.

SourceFiles
    volumeWeighted.C

\*---------------------------------------------------------------------------*/

#ifndef volumeWeighted_H
#define volumeWeighted_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volumeWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volumeWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const volumeWeighted&);


public:

    //- Runtime type information
    TypeName("volumeWeighted");


    // Constructors

        //- Construct from mesh
        volumeWeighted(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        volumeWeighted(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        volumeWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& sf
        ) const
        {
            const fvMesh& mesh = sf.mesh();
            const scalarField& V = mesh.V();
            const labelList& own = mesh.owner();
            const labelList& nei = mesh.neighbour();
        
            tmp<surfaceScalarField> taw
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "volumeWeightedWeights",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar("0.5", dimless, 0.5)
                )
            );
            // Modify the weights so that they are volume weighted
            surfaceScalarField& aw = taw();
            forAll(aw, faceI)
            {
                aw[faceI] = V[own[faceI]]/(V[own[faceI]] + V[nei[faceI]]);
            }

            surfaceScalarField::GeometricBoundaryField& awbf =
                taw().boundaryField();

            forAll(awbf, patchi)
            {
                if (!awbf[patchi].coupled())
                {
                    awbf[patchi] = 1.0;
                }
            }

            return taw;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
