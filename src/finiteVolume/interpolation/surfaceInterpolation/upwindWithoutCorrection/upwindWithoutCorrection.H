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
    Foam::upwindWithoutCorrection

Description
    upwindWithoutCorrection interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction. Third order accuracy on a regular grid. Otherwise 1st order. 
    Scheme is upwindWithoutCorrection below Courant number Co1 and upwind above Co2 and
    blended in between. If Co2 < Co1 or Co1 < 0 then no blending is done.

SourceFiles
    upwindWithoutCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef upwindWithoutCorrection_H
#define upwindWithoutCorrection_H

#include "upwind.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class upwindWithoutCorrection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindWithoutCorrection
:
    public upwind<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const upwindWithoutCorrection&) = delete;

public:

    //- Runtime type information
    TypeName("upwindWithoutCorrection");


    // Constructors

        //- Construct from faceFlux and Istream
        upwindWithoutCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        );

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        upwindWithoutCorrection
        (
            const fvMesh& mesh,
            Istream& is
        );

        //- Construct from faceFlux and Istream
        upwindWithoutCorrection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        upwindWithoutCorrection(const upwindWithoutCorrection&) = delete;


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


// Specialise for volVectorField to support cached gradients
template<>
tmp<surfaceVectorField> upwindWithoutCorrection<vector>::correction
(
    const volVectorField& vf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "upwindWithoutCorrectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
