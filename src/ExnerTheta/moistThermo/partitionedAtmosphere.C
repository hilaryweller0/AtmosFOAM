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
#include "surfaceInterpolate.H"
#include "fvcDiv.H"

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
        IOobject("rho", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("rho", dimDensity, scalar(0))
    ),
    dRhodt_
    (
        IOobject("dRhodt", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("dRhodt", dimDensity/dimTime, scalar(0))
    ),
    theta_
    (
        IOobject("theta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("theta", dimTemperature, scalar(0))
    ),
    Uf_
    (
        IOobject("Uf", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedVector("Uf", dimVelocity, vector::zero)
    ),
    flux_
    (
        IOobject("flux", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0))
    ),
    Psi_
    (
        IOobject("Psi", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("Psi", dimDensity, scalar(0))
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
    
    rho_.oldTime();
    dRhodt_.oldTime();
    flux_.oldTime();
    rho_.write();
    theta_.write();
    Uf_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::~partitionedAtmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::partitionedAtmosphere::sumDensity()
{
    // The zero'th partition
    const partition& part = operator[](0);

    // Initialise the density as density for partion 0
    rho_ = part.sigmaRho();

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        rho_ += part.sigmaRho();
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
    const partition& part = operator[](0);

    // Initialise the momentum from partion 0
    surfaceScalarField rhof = linearInterpolate(part.sigmaRho());
    surfaceVectorField rhoU = rhof*part.Uf();
    surfaceScalarField rhofSum = rhof;

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        const partition& part = operator[](ipart);
        rhof = linearInterpolate(part.sigmaRho());
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
    const perfectGasPhase& air = operator[](0).operator[](0).gas();
    const scalar kappa = air.kappa();

    // The averaging of theta over partitions is a complicated sum derived
    // so that the equation of state holds in each partition

    // The zero'th partition
    partition& part = operator[](0);

    // Initialise the sums from partion 0
    volScalarField thetaSum = part.sigma()*pow
    (
        part.rhoR()*part.theta()/(air.p0()*part.volGas()),
        kappa/(1-kappa)
    );

    // Sum contributions from other partitions
    for (label ipart = 1; ipart < size(); ipart++)
    {
        partition& part = operator[](ipart);
        thetaSum += part.sigma()*pow
        (
            part.rhoR()*part.theta()/(air.p0()*part.volGas()),
            kappa/(1-kappa)
        );
    }
    
    theta_ = air.p0()*(1-volLiquid())/rhoR()*pow(thetaSum, (1-kappa)/kappa);
    
    return theta_;
}

void Foam::partitionedAtmosphere::updateSigmas(const volScalarField& Exner)
{
    //const volScalarField ExnerForRho = sumDensity()/Psi();

    volScalarField sumSigma
    (
        IOobject("sumSigma", Exner.mesh().time().timeName(), Exner.mesh()),
        Exner.mesh(), dimensionedScalar("sum", dimless, scalar(0))
    );

    // Re-calculate all the sigmas
    for(label ip = 0; ip < size(); ip++)
    {
        partition& parti = operator[](ip);
        parti.updateSigma(Exner); //ForRho);
        sumSigma += parti.sigma();
    }
    
    Info << "1-sumSigma goes from " << 1-max(sumSigma).value() << " to "
         << 1-min(sumSigma).value() << endl;
         
    // Scale so that sigmas sum to 1
    for(label ip = 0; ip < size(); ip++)
    {
        partition& parti = operator[](ip);
        parti.sigma() /= sumSigma;
    }
}


Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::rhoR() const
{
    tmp<volScalarField> trhoRt
    (
        new volScalarField
        (
            IOobject
            (
                "rhoRt",
                rho_.time().timeName(),
                rho_.mesh()
            ),
            operator[](0).sigma()*operator[](0).rhoR()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoRt.ref() += operator[](ip).sigma()*operator[](ip).rhoR();
    }
    return trhoRt;
}


Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::volLiquid()
{
    operator[](0).sumVolLiquid();

    tmp<volScalarField> tvol
    (
        new volScalarField
        (
            IOobject
            (
                "volGas",
                rho_.time().timeName(),
                rho_.mesh()
            ),
            operator[](0).sigmaVolLiquid()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        operator[](ip).sumVolLiquid();
        tvol.ref() += operator[](ip).sigmaVolLiquid();
    }
    return tvol;
}


Foam::tmp<Foam::volScalarField> Foam::partitionedAtmosphere::ExnerFromState()
{
    const perfectGasPhase& air = operator[](0).operator[](0).gas();
    const scalar kappa = air.kappa();
    const fvMesh& mesh = air.rho().mesh();
    const Time& time = air.rho().time();

    tmp<volScalarField> tExner
    (
        new volScalarField
        (
            IOobject("Exner", time.timeName(), mesh),
            pow
            (
                rhoR()*updateTheta()/(air.p0()*(1-volLiquid())), 
                kappa/(1-kappa)
            )
        )
    );

    return tExner;
}


Foam::volScalarField& Foam::partitionedAtmosphere::ExnerFromState
(
    volScalarField& Exner
)
{
    const perfectGasPhase& air = operator[](0).operator[](0).gas();
    const scalar kappa = air.kappa();

    Exner == pow
    (
        rhoR()*updateTheta()/(air.p0()*(1-volLiquid())),
        kappa/(1-kappa)
    );

    return Exner;
}


void Foam::partitionedAtmosphere::setGradPcoeff
(
    surfaceScalarField& gradPcoeff
)
{
    const perfectGasPhase& air = operator[](0).operator[](0).gas();

    gradPcoeff == dimensionedScalar
    (
        "gradPcoeff",
        dimensionSet(1,-1,-2,0,0), scalar(0)
    );
    
    for (label ip = 0; ip < size(); ip++)
    {
        partition& parti = operator[](ip);
        gradPcoeff += air.Cp()*parti.updateThetaRho()
                      *fvc::interpolate(parti.sigmaRho());
    }
}


void Foam::partitionedAtmosphere::updateCompressibility
(
    const volScalarField& Exner
)
{
    sumDensity();
    dRhodt() = -fvc::div(flux());
    Psi_ == rho()/Exner;
    for (label ip = 0; ip < size(); ip++)
    {
        partition& parti = operator[](ip);
        parti.Psi() = parti.sigmaRho()/(parti.sigma()*Exner);
    }
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
    flux_.write();
    Uf_.write();
    static_cast<volScalarField>(ExnerFromState()).write();
}


void Foam::partitionedAtmosphere::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.readUpdate();
    }
    sumDensity();
    updateUf();
    updateFlux();
    updateTheta();
}


// ************************************************************************* //
