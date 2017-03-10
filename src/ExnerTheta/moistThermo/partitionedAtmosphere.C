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
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::partitionedAtmosphere
(
    const wordList& partitionNames,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    PtrList<partition>(partitionNames.size()),
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("rho", dimDensity, scalar(0))
    ),
    dRhodt_
    (
        IOobject("dRhodt", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("dRhodt", dimDensity-dimTime, scalar(0))
    ),
    theta_
    (
        IOobject("theta", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("theta", dimTemperature, scalar(0))
    ),
    Uf_
    (
        IOobject("Uf", mesh.time().timeName(), mesh),
        mesh, dimensionedVector("Uf", dimVelocity, vector::zero)
    ),
    flux_
    (
        IOobject("flux", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0))
    )
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
    sumDensity();
    updateUf();
    updateFlux();
    updateTheta();
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
    PtrList<partition>(partitionNames.size()),
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("rho", dimDensity, scalar(0))
    ),
    dRhodt_
    (
        IOobject("dRhodt", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("dRhodt", dimDensity/dimTime, scalar(0))
    ),
    theta_
    (
        IOobject("theta", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("theta", dimTemperature, scalar(0))
    ),
    Uf_
    (
        IOobject("Uf", mesh.time().timeName(), mesh),
        mesh, dimensionedVector("Uf", dimVelocity, vector::zero)
    ),
    flux_
    (
        IOobject("flux", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0))
    )
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
    
    sumDensity();
    updateUf();
    updateFlux();
    updateTheta();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::~partitionedAtmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::partitionedAtmosphere::sumDensity()
{
    // The zero'th partition
    partition& part = operator[](0);

    // Initialise the density as density for partion 0
    rho_ = part.sumDensity();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        partition& part = operator[](ipart);
        rho_ += part.sumDensity();
    }
    
    return rho_;
}

//Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::thetae() const
//{
//    // The total density
//    volScalarField rhoSum = rho();

//    // The zero'th partition
//    const partition& part = operator[](0);

//    // Initialise the thetae as thetae for partion 0
//    tmp<volScalarField> tt
//    (
//        new volScalarField
//        (
//            IOobject
//            (
//                "thetae",
//                part.T().time().timeName(),
//                part.T().mesh()
//            ),
//            part.thetae(part.T())*part.sumDensity()*part.sigma()/rhoSum
//        )
//    );
//    volScalarField& thetae = tt.ref();

//    // Sum contributions from other partitions
//    for (label ipart = 1; ipart < size(); ipart++)
//    {
//        const partition& part = operator[](ipart);
//        thetae += part.thetae(part.T())*part.sumDensity()*part.sigma()/rhoSum;
//    }
//    
//    return tt;
//}

Foam::surfaceVectorField& Foam::partitionedAtmosphere::updateUf()
{
    // The zero'th partition
    partition& part = operator[](0);

    // Initialise the momentum from partion 0
    surfaceScalarField rhof = linearInterpolate(part.sumDensity());
    surfaceVectorField rhoU = rhof*part.Uf();
    surfaceScalarField rhofSum = rhof;

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        partition& part = operator[](ipart);
        rhof = linearInterpolate(part.sumDensity());
        rhoU += rhof*part.Uf();
        rhofSum += rhof;
    }
    
    Uf_ = rhoU/rhofSum;
    return Uf_;
}


Foam::surfaceScalarField& Foam::partitionedAtmosphere::updateFlux()
{
    // The zero'th partition
    partition& part = operator[](0);

    // Initialise the momentum from partion 0
    flux_ = part.flux();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        partition& part = operator[](ipart);
        flux_ += part.flux();
    }
    
    return flux_;
}


Foam::volScalarField& Foam::partitionedAtmosphere::updateTheta()
{
    // The zero'th partition
    partition& part = operator[](0);

    // Initialise the momentum from partion 0
    volScalarField rho = part.sumDensity();
    volScalarField rhoTheta = rho*part.theta();
    volScalarField rhoSum = rho;

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        partition& part = operator[](ipart);
        rho = part.sumDensity();
        rhoTheta += rho*part.theta();
        rhoSum += rho;
    }
    
    theta_ = rhoTheta/rhoSum;
    return theta_;
}


void Foam::partitionedAtmosphere::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.write();
    }
    rho_.write();
    theta_.write();
}


void Foam::partitionedAtmosphere::readUpdate(const volScalarField& Exner)
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.readUpdate(Exner);
    }
    sumDensity();
}


// ************************************************************************* //
