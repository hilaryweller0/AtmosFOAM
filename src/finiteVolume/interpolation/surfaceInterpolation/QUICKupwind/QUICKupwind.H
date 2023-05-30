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
    Foam::QUICKupwind

Description
    QUICKupwind interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction. Same weighting factors as QUICK but all applied as an explicit
    correction on first-order upwind. 
    Scheme is QUICK below Courant number Co1 and upwind above Co2 and
    blended in between.

SourceFiles
    QUICKupwind.C

\*---------------------------------------------------------------------------*/

#ifndef QUICKupwind_H
#define QUICKupwind_H

#include "upwind.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class QUICKupwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class QUICKupwind
:
    public upwind<Type>
{
    // Private Data

        //- Courant number below which QUICK is used
        const scalar Co1_;

        //- Courant number above which upwind
        const scalar Co2_;

public:

    //- Runtime type information
    TypeName("QUICKupwind");


    // Constructors

        //- Construct from faceFlux and Istream
        QUICKupwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        );

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        QUICKupwind
        (
            const fvMesh& mesh,
            Istream& is
        );

        //- Construct from faceFlux and Istream
        QUICKupwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        QUICKupwind(const QUICKupwind&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Blending factor to scale the correction on upwind
        tmp<surfaceScalarField> blendingFactor() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const QUICKupwind&) = delete;
};


// Specialise for volVectorField to support cached gradients
template<>
tmp<surfaceVectorField> QUICKupwind<vector>::correction
(
    const volVectorField& vf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QUICKupwindTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
