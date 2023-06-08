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

#include "deformationalNonDivergentGeodesicVelocityFieldData.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformationalNonDivergentGeodesicVelocityFieldData, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::deformationalNonDivergentGeodesicVelocityFieldData
::deformationalNonDivergentGeodesicVelocityFieldData
(
    const fvMesh& mesh,
    const scalar radius,
    const scalar deformScale,
    const scalar endTime
)
:
    MeshObject
    <
        fvMesh, Foam::MoveableMeshObject,
        deformationalNonDivergentGeodesicVelocityFieldData
    >(mesh),
    radius_(radius),
    deformScale_(deformScale),
    endTime_(endTime),
    phat_sqrrByT_(mesh.nPoints()),
    Dcos2Lat_(mesh.nPoints()),
    TwoPiSinLat_(mesh.nPoints())
{
    calcData();
}


void Foam::deformationalNonDivergentGeodesicVelocityFieldData::calcData()
{
    scalarField magp = mag(mesh().points());

    phat_sqrrByT_ = mesh().points()/magp*sqr(radius_)/endTime_;

    Dcos2Lat_ = deformScale_*
    (
        1 - sqr(mesh().points().component(2)/magp)
    );

    TwoPiSinLat_ = 2*M_PI*mesh().points().component(2)/magp;
}

bool Foam::deformationalNonDivergentGeodesicVelocityFieldData::movePoints()
{
    calcData();
    return true;
}


// ************************************************************************* //
