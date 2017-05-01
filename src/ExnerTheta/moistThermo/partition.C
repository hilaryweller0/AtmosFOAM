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
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partition::partition
(
    const word& partitionName__,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    baseAtmosphere(add(partitionName__, partNames), mesh, dict),
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
    sigmaRho_
    (
        IOobject
        (
            partitionName_+"sigmaRho", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        sigma_*baseAtmosphere::sumDensity()
    ),
    theta_
    (
        IOobject
        (
            partitionName_+"theta", mesh.time().timeName(), mesh,
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
    Uf_
    (
        IOobject
        (
            partitionName_+"Uf", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    u_
    (
        IOobject(partitionName_+"u", mesh.time().timeName(), mesh),
        fvc::reconstruct(Uf_ & mesh.Sf()),
        "slip"
    ),
    flux_
    (
        IOobject
        (
            partitionName_+"flux", mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        linearInterpolate(sigmaRho_)*(Uf_ & mesh.Sf()),
        "fixedValue"
    ),
    dFluxdt_
    (
        IOobject(partitionName_+"dFluxdt", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("dFluxdt", dimensionSet(1,0,-2,0,0), scalar(0)),
        "fixedValue"
    ),
    sigmaVolLiquid_
    (
        IOobject(partitionName_+"sigmaVolLiquid", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("sigmaVolLiquid", dimless, scalar(0))
    ),
    dSigmaRhodt_
    (
        IOobject(partitionName_+"dSigmaRhodt", mesh.time().timeName(), mesh),
        -fvc::div(flux_)
    ),
    dSigmaRhoThetadt_
    (
        IOobject(partitionName_+"dSigmaRhoThetadt", mesh.time().timeName(), mesh),
        -fvc::div(flux_, theta_)
    ),
    thetaRho_
    (
        IOobject(partitionName_+"thetaRho", mesh.time().timeName(), mesh),
        fvc::interpolate(theta_)
    ),
    Psi_
    (
        IOobject(partitionName_+"Psi", mesh.time().timeName(), mesh),
        sigmaRho_
    )
{
    sumVolLiquid();
    T_ = theta_*ExnerFromState();
    updateThetaRho();
    sigmaRho_.write();

    sigma_.oldTime();
    sigmaRho_.oldTime();
    dSigmaRhodt_.oldTime();
    theta_.oldTime();
    Uf_.oldTime();
    flux_.oldTime();
    dFluxdt_.oldTime();
    dSigmaRhoThetadt_.oldTime();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partition::~partition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::partition::sumDensity()
{
    sigmaRho_ == dimensionedScalar("zero", dimDensity, scalar(0));
    for(label is = 0; is < size(); is++)
    {
        sigmaRho_ += operator[](is).gas().rho()
                  + operator[](is).liquid().rho()*operator[](is).liquid().v();
    }
    
    sigmaRho_ *= sigma_;
    
    return sigmaRho_;
}


Foam::volScalarField& Foam::partition::sumVolLiquid()
{
    sigmaVolLiquid_ = dimensionedScalar("zero", dimless, scalar(0));
    for(label is = 0; is < size(); is++)
    {
        sigmaVolLiquid_ += operator[](is).liquid().v();
    }
    sigmaVolLiquid_*= sigma_;
    return sigmaVolLiquid_;
}


Foam::volScalarField& Foam::partition::updateSigma(const volScalarField& Exner)
{
    sigma_ = sigmaRho_/(Psi_*Exner);
    return sigma_;
}


Foam::tmp<Foam::volScalarField> Foam::partition::ExnerFromState() const
{
    const perfectGasPhase& air = operator[](0).gas();
    const fvMesh& mesh = air.rho().mesh();
    const Time& time = air.rho().time();

    tmp<volScalarField> tExner
    (
        new volScalarField
        (
            IOobject("Exner", time.timeName(), mesh),
            pow
            (
                rhoR()*theta()/(air.p0()*volGas()), 
                air.kappa()/(1-air.kappa())
            )
        )
    );

    return tExner;
}


Foam::volScalarField& Foam::partition::ExnerFromState
(
    volScalarField& Exner
) const
{
    const perfectGasPhase& air = operator[](0).gas();

    Exner = pow
    (
        rhoR()*theta()/(air.p0()*volGas()), 
        air.kappa()/(1-air.kappa())
    );

    return Exner;
}


Foam::tmp<Foam::volScalarField> Foam::partition::thetaSource() const
{
    const perfectGasPhase& air = operator[](0).gas();
    const fvMesh& mesh = air.rho().mesh();
    const Time& time = air.rho().time();

    // Initialise the source term as -divu
    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "Stheta",
                time.timeName(),
                mesh
            ),
            -fvc::div(flux()/fvc::interpolate(sigmaRho()))
        )
    );
    volScalarField& S = tS.ref();
    
    // Pre-calculate rhoRt, rhoCp and rhoCv
    const volScalarField rhoR_ = rhoR();
    const volScalarField rhoCp_ = rhoCp();
    const volScalarField rhoCv_ = rhoCv();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Scale the divergence term
    S *= rhoR_/rhoCv_ - air.kappa()*rhoCp_/rhoCv_;
    
    // Add the terms relating to condensation for each species
    for(label ip = 0; ip < size(); ip++)
    {
        const fluidSpecie& species = operator[](ip);
    
        // Only if there is any condensation
        if (species.pvs0() < species.gas().p0())
        {
            S += species.condensation()/dt/rhoCv_*
             (
                 air.Cv()*species.latentHeat(T())/(air.Cp()*T())
               - species.gas().R()*(1 - air.kappa()*rhoCp_/rhoR_)
             );
         }
    }
    
    return tS;
}

Foam::surfaceScalarField& Foam::partition::updateThetaRho()
{
    const perfectGasPhase& air = operator[](0).gas();
    thetaRho_ = fvc::interpolate
    (
        sigma()*rhoR()*theta()/(air.R()*sigmaRho()*volGas())
    );
    return thetaRho_;
}


void Foam::partition::write()
{
    baseAtmosphere::write();
    sigma_.write();
    sigmaRho_.write();
    T_.write();
    theta_.write();
    Uf_.write();
    flux_.write();
}


void Foam::partition::readUpdate()
{
    const fvMesh& mesh = sigma_.mesh();
    baseAtmosphere::readUpdate();
    sigma_ = volScalarField
    (
        IOobject(partitionName_+"sigma", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        sigma_
    );
    theta_ = volScalarField
    (
        IOobject(partitionName_+"theta", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        theta_
    );
    Uf_ = surfaceVectorField
    (
        IOobject(partitionName_+"Uf", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        Uf_
    );
    u_ = fvc::reconstruct(Uf_ & mesh.Sf());
    flux_ = linearInterpolate(sumDensity())*(Uf_ & mesh.Sf());
    sumDensity();
    sumVolLiquid();
    T_ = theta_*ExnerFromState();
}

// ************************************************************************* //
