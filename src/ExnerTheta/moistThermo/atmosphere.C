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

#include "atmosphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphere::atmosphere
(
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    PtrList<fluidSpecie>(partNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new fluidSpecie
            (
                partNames[ip],
                IOobject(partNames[ip]+"VapourRho", mesh.time().timeName(), mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                IOobject(partNames[ip]+"LiquidFrac", mesh.time().timeName(), mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                mesh,
                dict.subDict(partNames[ip])
            )
        );
    }
    sumDensity();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphere::~atmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::atmosphere::volGas() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tvG
    (
        new volScalarField
        (
            IOobject
            (
                "volGas",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            1 - operator[](0).liquid().v()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        tvG.ref() -= operator[](ip).liquid().v();
    }
    return tvG;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::volAir() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tvA
    (
        new volScalarField
        (
            IOobject
            (
                "volAir",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            volGas()*air.R()/rhoR()*air.rho()
        )
    );
    return tvA;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::sumDensity() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tRho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
        )
    );

    for(label ip = 1; ip < size(); ip++)
    {
        tRho.ref() += operator[](ip).gas().rho()
                   + operator[](ip).liquid().rho()*operator[](ip).liquid().v();
    }
    return tRho;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoR() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoRt
    (
        new volScalarField
        (
            IOobject
            (
                "rhoRt",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.R()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoRt.ref() += operator[](ip).gas().rho()*operator[](ip).gas().R();
    }
    return trhoRt;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoCp() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoCp
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCp",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.Cp()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
            *operator[](0).liquid().Cp()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoCp.ref() += operator[](ip).gas().rho()*operator[](ip).gas().Cp()
                    + operator[](ip).liquid().rho()*operator[](ip).liquid().v()
                      *operator[](ip).liquid().Cp();
    }
    return trhoCp;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoCv() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoCv
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCv",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.Cv()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
            *operator[](0).liquid().Cp()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoCv.ref() += operator[](ip).gas().rho()*operator[](ip).gas().Cv()
                    + operator[](ip).liquid().rho()*operator[](ip).liquid().v()
                      *operator[](ip).liquid().Cp();
    }
    return trhoCv;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::sumPressure
(
    const volScalarField& T
) const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tp
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.partialPressure(T)
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        tp.ref() += operator[](ip).gas().partialPressure(T);
    }
    // Divide by the volume fraction occupied by gas
    tp.ref() /= volGas();
    return tp;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::pFromT
(
    const volScalarField& T
) const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tp
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            rhoR()*T/volGas()
        )
    );
    return tp;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::ExnerFromTheta
(
    const volScalarField& theta
) const
{
    const perfectGasPhase& air = operator[](0).gas();
    const scalar kappa = air.kappa();
    const dimensionedScalar& p0 = air.p0();
    
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "Exner",
                theta.time().timeName(),
                theta.mesh()
            ),
            pow(theta*rhoR()/(p0*volGas()), kappa/(1-kappa)),
            wordList(theta.boundaryField().size(), "hydrostaticExner")
        )
    );
    return tE;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoFromP
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                p.time().timeName(),
                p.mesh()
            ),
            p*volGas()/(rhoR()/sumDensity()*T)
        )
    );
    return trho;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoFromExner
(
    const volScalarField& Exner,
    const volScalarField& theta
) const
{
    const perfectGasPhase& air = operator[](0).gas();
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                Exner.time().timeName(),
                Exner.mesh()
            ),
            pow(Exner, (1-air.kappa())/air.kappa())*air.p0()*volGas()
           /(theta*rhoR()/sumDensity())
        )
    );
    return trho;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::thetaSource
(
    const volScalarField& T,
    const volScalarField& divu
) const
{
    const perfectGasPhase& air = operator[](0).gas();

    // Initialise the source term as -divu
    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "Stheta",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            -divu
        )
    );
    volScalarField& S = tS.ref();
    
    // Pre-calculate rhoRt, rhoCp and rhoCv
    const volScalarField rhoR_ = rhoR();
    const volScalarField rhoCp_ = rhoCp();
    const volScalarField rhoCv_ = rhoCv();
    const dimensionedScalar& dt = S.mesh().time().deltaT();
    
    // Scale the divergence term
    S *= rhoR_/rhoCv_ - air.kappa()*rhoCp_/rhoCv_;
    
    // Add the terms relating to condensation for each phase
    for(label ip = 0; ip < size(); ip++)
    {
        const fluidSpecie& phase = operator[](ip);
    
        // Only if there is any condensation
        if (phase.pvs0() < phase.gas().p0())
        {
            S += phase.condensation()/dt/rhoCv_*
             (
                 air.Cv()*phase.latentHeat(T)/(air.Cp()*T)
               - phase.gas().R()*(1 - air.kappa()*rhoCp_/rhoR_)
             );
         }
    }
    
    return tS;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::thetae
(
    const volScalarField& T
) const
{
    const perfectGasPhase& air = operator[](0).gas();
    const fluidSpecie& water = operator[](1);

    // Initialise the thetae as T
    tmp<volScalarField> tt
    (
        new volScalarField
        (
            IOobject
            (
                "thetae",
                T.time().timeName(),
                T.mesh()
            ),
            T
        )
    );
    volScalarField& thetae = tt.ref();
    
    volScalarField Cp = air.Cp()
                      + water.liquid().Cp()*(sumDensity()/air.rho()-1);
    volScalarField a = pow(air.partialPressure(T)/air.p0(), -air.R()/Cp);
    volScalarField b = water.gas().rho()*water.latentHeat(T)/(Cp*air.rho());

    // The residual to miniminse every Newton step
    thetae = T*a*exp(b/T);
    
    return tt;
}

void Foam::atmosphere::TfromThetae
(
    volScalarField& T,
    const dimensionedScalar& thetae0,
    const scalar rt
) const
{
    const perfectGasPhase& air = operator[](0).gas();
    const fluidSpecie& water = operator[](1);

    // Setup variables for a Newton method
    dimensionedScalar Cp = air.Cp() + water.liquid().Cp()*rt;
    volScalarField a = pow(air.partialPressure(T)/air.p0(), -air.R()/Cp);
    volScalarField b = water.gas().rho()*water.latentHeat(T)/(Cp*air.rho());

    // The residual to miniminse every Newton step
    volScalarField resid = thetae0 - T*a*exp(b/T);
    scalar RMSresid = Foam::sqrt(sum(sqr(resid.internalField())).value());
    
    // convergence criterial
    scalar conv = 1e-6;
    
    // Newton steps
    for(label i = 0; i < 10 && RMSresid > conv; i++)
    {
        T -= resid*T/(a*(b-T)*exp(b/T));
        resid = thetae0 - T*a*exp(b/T);
        RMSresid = Foam::sqrt(sum(sqr(resid.internalField())).value());
    }
}

void Foam::atmosphere::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        const fluidSpecie& phase = operator[](ip);
        phase.gas().rho().write();
        phase.liquid().v().write();
    }
}


void Foam::atmosphere::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        fluidSpecie& phase = operator[](ip);
        const fvMesh& mesh = phase.gas().rho().mesh();
        
        phase.gas().rho() = volScalarField
        (
            IOobject(phase.name()+"VapourRho", mesh.time().timeName(), mesh,
                     IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            phase.gas().rho()
        );
        phase.liquid().v() = volScalarField
        (
            IOobject(phase.name()+"LiquidFrac", mesh.time().timeName(), mesh,
                     IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            phase.liquid().v()
        );
    }
}


// ************************************************************************* //
