/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::quinticUpwind

Description
    quinticUpwind interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction. Fifth order accuracy on a regular grid. Otherwise 1st order.

SourceFiles
    quinticUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef quinticUpwind_H
#define quinticUpwind_H

#include "upwind.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class quinticUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quinticUpwind
:
    public upwind<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const quinticUpwind&) = delete;

public:

    //- Runtime type information
    TypeName("quinticUpwind");


    // Constructors

        //- Construct from faceFlux and Istream
        quinticUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        );

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        quinticUpwind
        (
            const fvMesh& mesh,
            Istream& is
        );

        //- Construct from faceFlux and Istream
        quinticUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        quinticUpwind(const quinticUpwind&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>&
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "quinticUpwindTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
