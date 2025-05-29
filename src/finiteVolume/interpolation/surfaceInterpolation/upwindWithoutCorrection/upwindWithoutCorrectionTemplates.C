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

\*---------------------------------------------------------------------------*/

#include "upwindWithoutCorrection.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::upwindWithoutCorrection<Type>::upwindWithoutCorrection
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux
)
:
    upwind<Type>(mesh, faceFlux)
{}


template<class Type>
Foam::upwindWithoutCorrection<Type>::upwindWithoutCorrection
(
    const fvMesh& mesh,
    Istream& is
)
:
    upwind<Type>(mesh, is)
{}


template<class Type>
Foam::upwindWithoutCorrection<Type>::upwindWithoutCorrection
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    upwind<Type>(mesh, faceFlux, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::upwindWithoutCorrection<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "upwindWithoutCorrection::correction(" + vf.name() + ')',
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    return tsfCorr;
}

// ************************************************************************* //
