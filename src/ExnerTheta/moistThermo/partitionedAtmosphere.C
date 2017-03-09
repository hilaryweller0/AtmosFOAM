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

#include "partitionedAtmosphere.H"
#include "moreListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::partitionedAtmosphere
(
    const wordList& partitionNames,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    PtrList<partition>(partitionNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new partition
            (
                partitionNames[ip],
                partNames,
                mesh,
                dict
            )
        );
    }
}


Foam::partitionedAtmosphere::partitionedAtmosphere
(
    const wordList& partitionNames,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict,
    const volScalarField& Exner
)
:
    PtrList<partition>(partitionNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new partition
            (
                partitionNames[ip],
                partNames,
                mesh,
                dict,
                Exner
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::~partitionedAtmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::rho() const
{
    // The zero'th partition
    const partition& part = operator[](0);

    // Initialise the density as density for partion 0
    tmp<volScalarField> rhot
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                part.T().time().timeName(),
                part.T().mesh()
            ),
            part.sumDensity()*part.sigma()
        )
    );
    volScalarField& rhoSum = rhot.ref();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        rhoSum += part.sumDensity()*part.sigma();
    }
    
    return rhot;
}

Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::thetae() const
{
    // The total density
    volScalarField rhoSum = rho();

    // The zero'th partition
    const partition& part = operator[](0);

    // Initialise the thetae as thetae for partion 0
    tmp<volScalarField> tt
    (
        new volScalarField
        (
            IOobject
            (
                "thetae",
                part.T().time().timeName(),
                part.T().mesh()
            ),
            part.thetae(part.T())*part.sumDensity()*part.sigma()/rhoSum
        )
    );
    volScalarField& thetae = tt.ref();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        thetae += part.thetae(part.T())*part.sumDensity()*part.sigma()/rhoSum;
    }
    
    return tt;
}

Foam::tmp<Foam::surfaceVectorField> Foam::partitionedAtmosphere::rhoUf() const
{
    // The zero'th partition
    const partition& part = operator[](0);

    // Initialise the momentum from partion 0
    tmp<surfaceVectorField> rhoUt
    (
        new surfaceVectorField
        (
            IOobject
            (
                "rhoUf",
                part.T().time().timeName(),
                part.T().mesh()
            ),
            linearInterpolate(part.sumDensity()*part.sigma())*part.Uf()
        )
    );
    volScalarField& rhoU = rhoUt.ref();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        rhoU += linearInterpolate(part.sumDensity()*part.sigma())*part.Uf();
    }
    
    return rhoUt;
}

Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::rhoTheta() const
{
    // The zero'th partition
    const partition& part = operator[](0);

    // Initialise the momentum from partion 0
    tmp<volScalarField> rhoTt
    (
        new volScalarField
        (
            IOobject
            (
                "rhoTheta",
                part.T().time().timeName(),
                part.T().mesh()
            ),
            part.sumDensity()*part.sigma()*part.theta()
        )
    );
    volScalarField& rhoT = rhoTt.ref();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        rhoT += part.sumDensity()*part.sigma()*part.theta();
    }
    
    return rhoTt;
}

void Foam::partitionedAtmosphere::updateThetaT(const volScalarField& Exner)
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.updateThetaT(Exner);
    }
}


void Foam::partitionedAtmosphere::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.write();
    }
}


void Foam::partitionedAtmosphere::readUpdate(const volScalarField& Exner)
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.readUpdate(Exner);
    }
}


// ************************************************************************* //
