/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "liquidPhase.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidPhase::liquidPhase
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Cp_(dict.lookup("Cp")),
    rho_(dict.lookup("rho")),
    v_(io, mesh),
    dvdt_
    (
        IOobject("d"+v_.name()+"dt", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("", v_.dimensions()/dimTime, scalar(0))
    )
{
    v_.oldTime();
    dvdt_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidPhase::~liquidPhase()
{}

// ************************************************************************* //
