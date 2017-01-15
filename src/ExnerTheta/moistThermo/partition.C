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

#include "partition.H"
#include "moreListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partition::partition
(
    const word& partitionName__,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    atmosphere(add(partitionName__, partNames), mesh, dict),
    partitionName_(partitionName__),
    sigma_
    (
        IOobject
        (
            partitionName_+"sigma", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            partitionName_+"T", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, scalar(0))
    ),
    theta_
    (
        IOobject
        (
            partitionName_+"theta", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, scalar(0))
    )
{}


Foam::partition::partition
(
    const word& partitionName__,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict,
    const volScalarField& Exner
)
:
    atmosphere(add(partitionName__, partNames), mesh, dict),
    partitionName_(partitionName__),
    sigma_
    (
        IOobject
        (
            partitionName_+"sigma", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            partitionName_+"T", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, scalar(0))
    ),
    theta_
    (
        IOobject
        (
            partitionName_+"theta", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, scalar(0))
    )
{
    updateThetaT(Exner);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partition::~partition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::partition::updateThetaT(const volScalarField& Exner)
{
    // Update theta and T from Exner using the equation of state for the
    // perfect gas part of the partition
    perfectGasPhase& air = operator[](0).gas();
    const scalar kappa = air.kappa();
    const dimensionedScalar& p0 = air.p0();
    theta_ = p0*volGas()*pow(Exner, (1-kappa)/kappa)/rhoR();
    T_ = theta_*Exner;
}


void Foam::partition::write()
{
    atmosphere::write();
    sigma_.write();
    T_.write();
    theta_.write();
}


void Foam::partition::readUpdate(const volScalarField& Exner)
{
    atmosphere::readUpdate();
    updateThetaT(Exner);
}

// ************************************************************************* //
