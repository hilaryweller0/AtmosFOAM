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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::surfaceVectorField>
Foam::upwindWithoutCorrection<Foam::vector>::correction
(
    const volVectorField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceVectorField> tsfCorr
    (
        surfaceVectorField::New
        (
            "upwindWithoutCorrection::correction(" + vf.name() + ')',
            mesh,
            dimensioned<vector>(vf.name(), vf.dimensions(), Zero)
        )
    );

    return tsfCorr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makelimitedSurfaceInterpolationScheme(upwindWithoutCorrection)
}

// ************************************************************************* //
