/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    totalVariation

Description
    Calculates and the total variation of a volScalarField.

\*---------------------------------------------------------------------------*/

#ifndef totalVariation_H
#define totalVariation_H

namespace Foam
{
inline scalar totalVariation
(
    const volScalarField& vf
)
{
    const fvMesh& mesh = vf.mesh();
    scalar TV = 0;

    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];
        TV += mag(vf[own] - vf[nei]);
    }

    return TV;
}
}
#endif

