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

#include "fluidSpecie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSpecie::fluidSpecie
(
    const word& name__,
    const IOobject& gasIO,
    const IOobject& liquidIO,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name__),
    gas_(gasIO, mesh, dict),
    liquid_(liquidIO, mesh, dict.subDict("liquidDict")),
    Lv0_(dict.lookup("Lv0")),
    pvs0_(dict.lookup("pvs0")),
    condensation_("condensationOf"+name_, 0*gas_.rho())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSpecie::~fluidSpecie()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::fluidSpecie::updateCondensation
(
    const volScalarField& T
)
{
    volScalarField Sf = (gas_.partialPressure(T) - pSat(T))/(gas_.R()*T);
    volScalarField Sl = liquid_.v() * liquid_.rho();
    condensation_ = 0.5*(Sf - Sl + sqrt(sqr(Sf) + sqr(Sl)));
    return condensation_;
}


// ************************************************************************* //
