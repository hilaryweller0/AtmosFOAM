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

#include "perfectGasPhase.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::perfectGasPhase::perfectGasPhase
(
    const IOobject& io,
    const fvMesh& mesh,
    Istream& is,
    const dictionary& dict
)
:
    constTransport<hConstThermo<perfectGas<specie> > >(is),
    R_("R", dimGasConstant, constTransport<hConstThermo<perfectGas<specie> > >::R()),
    T0_(dict.lookup("T0")),
    p0_(dict.lookup("p0")),
    Cp_("Cp", dimGasConstant, cp(p0_.value(),T0_.value())/W()),
    Cv_(Cp_-R_),
    kappa_(R_.value()/Cp_.value()),
    rho_(io, mesh),
    dRhodt_
    (
        IOobject("d"+rho_.name()+"dt", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("", rho_.dimensions()/dimTime, scalar(0))
    )
{
    rho_.oldTime();
    dRhodt_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::perfectGasPhase::~perfectGasPhase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::perfectGasPhase::partialPressure
(
    const volScalarField& T
)
{
    tmp<volScalarField> tpp
    (
        new volScalarField
        (
            IOobject("pp", T.time().timeName(), T.mesh()),
            rho_*R()*T
        )
    );
    return tpp;
}

// ************************************************************************* //
