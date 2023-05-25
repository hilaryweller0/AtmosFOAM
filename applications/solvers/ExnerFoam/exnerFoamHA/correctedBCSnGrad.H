/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fv::correctedBCSnGrad

Description
    Simple central-difference snGrad scheme with non-orthogonal correction and
    second order extrapolation at external boundaries

SourceFiles
    correctedBCSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef correctedBCSnGrad_H
#define correctedBCSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class correctedBCSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correctedBCSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const correctedBCSnGrad&) = delete;


public:

    //- Runtime type information
    TypeName("correctedBC");


    // Constructors

        //- Construct from mesh
        correctedBCSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        correctedBCSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~correctedBCSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the correctedBCSnGrad
        //  for the given field using the gradient of the field
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        fullGradCorrection
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the explicit correction to the correctedBCSnGrad
        //  for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<surfaceScalarField> correctedBCSnGrad<scalar>::correction
(
    const volScalarField& vsf
) const;


template<>
tmp<surfaceVectorField> correctedBCSnGrad<vector>::correction
(
    const volVectorField& vvf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "correctedBCSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
