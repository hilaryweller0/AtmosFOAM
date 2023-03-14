/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "stratifiedRealizableKE.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> stratifiedRealizableKE<BasicTurbulenceModel>::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    const volScalarField W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
        )
    );

    tS.clear();

    const volScalarField phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    const volScalarField As(sqrt(6.0)*cos(phis));
    const volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    return 1.0/(A0_ + As*Us*k_/epsilon_);
}


template<class BasicMomentumTransportModel>
void stratifiedRealizableKE<BasicMomentumTransportModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    Cmu_ = rCmu(gradU, S2, magS);
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void stratifiedRealizableKE<BasicMomentumTransportModel>::correctNut()
{
    const volTensorField gradU(fvc::grad(this->U_));
    const volScalarField S2(modelName("S2"), 2*magSqr(dev(symm(gradU))));
    const volScalarField magS(modelName("magS"), sqrt(S2));

    correctNut(gradU, S2, magS);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
stratifiedRealizableKE<BasicMomentumTransportModel>::stratifiedRealizableKE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            1.46
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.44
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.08
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.2
        )
    ),
    sigmaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaTheta",
            this->coeffDict_,
            0.74
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    N2_
    (
        IOobject("N2", this->runTime_.timeName(),this->mesh_),
        this->mesh_,
        dimensionedScalar("N2", dimensionSet(0,0,-2,0,0), scalar(0))
    ),
    Ri_
    (
        IOobject("Ri", this->runTime_.timeName(),this->mesh_),
        this->mesh_,
        dimensionedScalar("Ri", dimless, scalar(0))
    ),
    Cmu_
    (
        IOobject("Cmu", this->runTime_.timeName(),this->mesh_),
        this->mesh_,
        dimensionedScalar("Cmu", dimless, scalar(0.09))
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    N2_.writeOpt() = IOobject::AUTO_WRITE;
    Ri_.writeOpt() = IOobject::AUTO_WRITE;
    Cmu_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool stratifiedRealizableKE<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        A0_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaTheta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void stratifiedRealizableKE<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const fvMesh& mesh = this->mesh_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const volScalarField& dbdz = mesh.objectRegistry::template
                 lookupObject<volScalarField>("dbdz");
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        modelName("divU"),
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(modelName("S2"), 2*magSqr(dev(symm(gradU))));
    const volScalarField magS(modelName("magS"), sqrt(S2));

    Cmu_ = rCmu(gradU, S2, magS);

    const volScalarField::Internal eta
    (
        modelName("eta"), magS()*k_()/epsilon_()
    );
    const volScalarField::Internal C1
    (
        modelName("C1"),
        max(eta/(scalar(5) + eta), scalar(0.43))
    );

    const volScalarField::Internal G
    (
        this->GName(),
        nut*(gradU.v() && dev(twoSymm(gradU.v())))
    );

    // Brunt-Vaisala frequency squared
    N2_ = dbdz;
    bound(N2_, dimensionedScalar("", dbdz.dimensions(), scalar(0)));

    // Set ghat (gravity direction)
    dimensionedVector ghat("ghat", dimless, vector(0,0,1));

    // Richardson number
    volScalarField magSqrdUdz = magSqr(ghat & gradU);
    Ri_ = N2_/bound
    (
        magSqrdUdz,
        dimensionedScalar("", dbdz.dimensions(), VSMALL)
    );
    
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
/*
    volScalarField epsDiffusion
    (
        "epsDiffusion",
        fvc::laplacian(alpha*rho*DepsilonEff(), epsilon_)/epsilon_
    );
    epsDiffusion.write();
*/
    surfaceScalarField snGradEpsilon
    (
        IOobject("snGradEpsilon", this->runTime_.timeName(),this->mesh_),
        fvc::snGrad(epsilon_)
    );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1*alpha()*rho()*magS()*epsilon_()
      - fvm::Sp
        (
            C2_*alpha()*rho()*epsilon_()/(k_() + sqrt(this->nu()()*epsilon_())),
            epsilon_
        )
      - C3_*Cmu_*N2_/sigmaTheta_*k_
      + C4_*min(sqrt(Ri_/C5_), scalar(1))*sqrt(N2_)*epsilon_
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    surfaceScalarField snGradk
    (
        IOobject("snGradk", this->runTime_.timeName(),this->mesh_),
        fvc::snGrad(k_)
    );

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      - fvm::Sp
        (
            Cmu_*N2_/sigmaTheta_*k_/epsilon_,
            k_
        )
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut(gradU, S2, magS);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
