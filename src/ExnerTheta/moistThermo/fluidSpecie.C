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
    gas_(gasIO, mesh, dict.lookup("gas"), dict.subDict("gasDict")),
    liquid_(liquidIO, mesh, dict.subDict("liquidDict")),
    Lv0_(dict.lookup("Lv0")),
    pvs0_(dict.lookup("pvs0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSpecie::~fluidSpecie()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidSpecie::latentHeat
(
    const volScalarField& T
)
{
    tmp<volScalarField> tLv
    (
        new volScalarField
        (
            IOobject("Lv", T.time().timeName(), T.mesh()),
            Lv0_ - (liquid_.Cp() - gas_.Cp())*(T - gas_.T0())
        )
    );
    return tLv;
}

Foam::tmp<Foam::volScalarField> Foam::fluidSpecie::pSat
(
    const volScalarField& T
)
{
    tmp<volScalarField> tpSat
    (
        new volScalarField
        (
            IOobject("pSat", T.time().timeName(), T.mesh()),
            pvs0_*exp(-Lv0_/gas_.R()*(1/T - 1/gas_.T0()))
        )
    );
    return tpSat;
}

Foam::tmp<Foam::volScalarField> Foam::fluidSpecie::condensation
(
    const volScalarField& T
)
{
    volScalarField Sf = (gas_.partialPressure(T) - pSat(T))/(gas_.R()*T);
    volScalarField Sl = liquid_.v() * liquid_.rho();

    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject("Scond", T.time().timeName(), T.mesh()),
            0.5*(Sf - Sl + sqrt(sqr(Sf) + sqr(Sl)))
        )
    );
    return tS;
}


// ************************************************************************* //
