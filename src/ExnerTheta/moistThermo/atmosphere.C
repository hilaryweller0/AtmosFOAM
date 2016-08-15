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
    PtrList<fluidSpecie>(partNames.size()),
    rho_
    (
        IOobject("rho", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("rho", dimensionSet(1,-3,0,0,0), scalar(0))
    )
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

Foam::tmp<Foam::volScalarField> Foam::atmosphere::volGas()
{
    perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tvG
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
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

Foam::volScalarField& Foam::atmosphere::sumDensity()
{
    perfectGasPhase& air = operator[](0).gas();

    rho_ = air.rho()
         + operator[](0).liquid().rho()*operator[](0).liquid().v();

    for(label ip = 1; ip < size(); ip++)
    {
        rho_ += operator[](ip).gas().rho()
              + operator[](ip).liquid().rho()*operator[](ip).liquid().v();
    }
    return rho_;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoRt()
{
    perfectGasPhase& air = operator[](0).gas();

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

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoCp()
{
    perfectGasPhase& air = operator[](0).gas();

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

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoCv()
{
    perfectGasPhase& air = operator[](0).gas();

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
)
{
    perfectGasPhase& air = operator[](0).gas();

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

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoThetaRho
(
    const volScalarField& T
)
{
    perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rhoThetaRho",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho() * sumPressure(T) / air.partialPressure(T)
        )
    );
    return trho;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::pFromT
(
    const volScalarField& T
)
{
    perfectGasPhase& air = operator[](0).gas();

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
            rhoRt()*T/volGas()
        )
    );
    return tp;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::ExnerFromTheta
(
    const volScalarField& theta
)
{
    perfectGasPhase& air = operator[](0).gas();
    const scalar kappa = air.kappa();
    const dimensionedScalar& p0 = air.p0();
    
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "Exner",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            pow(theta*rhoRt()/(p0*volGas()), kappa/(1-kappa)),
            wordList(theta.boundaryField().size(), "hydrostaticExner")
        )
    );
    return tE;
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::rhoFromP
(
    const volScalarField& p,
    const volScalarField& T
)
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
            p*volGas()/(rhoRt()/sumDensity()*T)
        )
    );
    return trho;
}

//Foam::tmp<Foam::volScalarField> Foam::atmosphere::ExnerFromThetaRho
//(
//    const volScalarField& ?
//)

Foam::tmp<Foam::volScalarField> Foam::atmosphere::thetaSource
(
    const volScalarField& T,
    const volScalarField& divu
)
{
    perfectGasPhase& air = operator[](0).gas();

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
    const volScalarField rhoRt_ = rhoRt();
    const volScalarField rhoCp_ = rhoCp();
    const volScalarField rhoCv_ = rhoCv();
    
    // Scale the divergence term
    S *= rhoRt_/rhoCv_ - air.R()/air.Cp()*rhoCp_/rhoCv_;
    
    // Add the terms relating to condensation for each phase
    for(label ip = 0; ip < size(); ip++)
    {
        fluidSpecie& phase = operator[](ip);
    
        // Only if there is any condensation
        if (phase.pvs0() < phase.gas().p0())
        {
            S += phase.condensation()/S.mesh().time().deltaT()/rhoCv_*
             (
                 air.Cv()*phase.latentHeat(T)/(air.Cp()*T)
               - phase.gas().R()*(1 - air.R()/air.Cp()*rhoCp_/rhoRt_)
             );
         }
    }
    
    return tS;
}

void Foam::atmosphere::write()
{
    rho().write();
    for(label ip = 0; ip < size(); ip++)
    {
        fluidSpecie& phase = operator[](ip);
        phase.gas().rho().write();
        phase.liquid().v().write();
    }
}

// ************************************************************************* //
