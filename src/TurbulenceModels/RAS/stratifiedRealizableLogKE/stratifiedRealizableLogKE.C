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

#include "stratifiedRealizableLogKE.H"
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
tmp<volScalarField> stratifiedRealizableLogKE<BasicTurbulenceModel>::rCmu
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
void stratifiedRealizableLogKE<BasicMomentumTransportModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    Cmu_ = rCmu(gradU, S2, magS);
    this->nut_ = Cmu_*sqr(k1)/eps1*exp(2*lnk_ - lnepsilon_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void stratifiedRealizableLogKE<BasicMomentumTransportModel>::correctNut()
{
    const volTensorField gradU(fvc::grad(this->U_));
    const volScalarField S2(typedName("S2"), 2*magSqr(dev(symm(gradU))));
    const volScalarField magS(typedName("magS"), sqrt(S2));

    correctNut(gradU, S2, magS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
stratifiedRealizableLogKE<BasicMomentumTransportModel>::stratifiedRealizableLogKE
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
    lnk_
    (
        IOobject
        (
            "lnk",
            this->runTime_.timeName(),this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    lnepsilon_
    (
        IOobject
        (
            "lnepsilon",
            this->runTime_.timeName(),this->mesh_,
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
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    Ri_.writeOpt() = IOobject::AUTO_WRITE;
    Cmu_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool stratifiedRealizableLogKE<BasicMomentumTransportModel>::read()
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
void stratifiedRealizableLogKE<BasicMomentumTransportModel>::correct()
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
        typedName("divU"),
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(typedName("S2"), 2*magSqr(dev(symm(gradU))));
    const volScalarField magS(typedName("magS"), sqrt(S2));

    Cmu_ = rCmu(gradU, S2, magS);

    const volScalarField::Internal eta
    (
        typedName("eta"), magS()*k_()/epsilon_()
    );
    const volScalarField::Internal C1
    (
        typedName("C1"),
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
    
    // Update lnepsilon, epsilon and G at the wall
    lnepsilon_.boundaryFieldRef().updateCoeffs();
    
    surfaceScalarField snGradlnepsilon
    (
        IOobject("snGradlnepsilon", this->runTime_.timeName(),this->mesh_),
        fvc::snGrad(lnepsilon_)
    );

/*    volScalarField epsDiffusion
    (
        "epsDiffusion",
        fvc::laplacian(DepsilonEff(), lnepsilon_)
      + min(DepsilonEff()*magSqr(fvc::grad(lnepsilon_)), 1/dt)
//      + DepsilonEff()*fvc::average(magSqr(fvc::snGrad(lnepsilon_)))
//      + DepsilonEff()*fvc::laplacian(lnepsilon_,lnepsilon_)
//      - DepsilonEff()*lnepsilon_*fvc::laplacian(lnepsilon_)
    );
    epsDiffusion.write();
    volScalarField epsSource
    (
        "epsSource",
        C1*alpha*rho*magS
      - C2_*Cmu_*alpha*rho*k_/(this->nut_ + sqrt(this->nut_*this->nu()()))
      - C3_*Cmu_*N2_/sigmaTheta_*k1/eps1*exp(lnk_-lnepsilon_)
      + C4_*min(sqrt(Ri_/C5_), scalar(1))*sqrt(N2_)
    );
    epsSource.write();
*/
    const dimensionedScalar& dt = this->runTime_.deltaT();
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, lnepsilon_)
      + fvm::div(alphaRhoPhi, lnepsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), lnepsilon_)
      - min(DepsilonEff()*magSqr(fvc::grad(lnepsilon_)), 1/dt)
      //- min(DepsilonEff()*fvc::average(magSqr(fvc::snGrad(lnepsilon_))), 1/dt)
      //- fvc::laplacian(alpha*rho*DepsilonEff(), epsilon_)/epsilon_
//      - alpha*rho*DepsilonEff()*fvc::laplacian(lnepsilon_, lnepsilon_)
//      + alpha*rho*DepsilonEff()*lnepsilon_*fvm::laplacian(lnepsilon_)
     ==
        C1*alpha()*rho()*magS()
      - C2_*Cmu_()*alpha()*rho()*k_()/(this->nut_() + sqrt(this->nut_()*this->nu()()))
      - C3_*Cmu_()*N2_()/sigmaTheta_*k1/eps1*exp(lnk_()-lnepsilon_())
      + C4_*min(sqrt(Ri_()/C5_), scalar(1))*sqrt(N2_())
      + fvModels.source(alpha, rho, lnepsilon_)
    );

    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(lnepsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    //bound(lnepsilon_, log(this->epsilonMin_/eps1));
    epsilon_ == eps1*exp(lnepsilon_);

    // Turbulent kinetic energy equation
    surfaceScalarField snGradlnk
    (
        IOobject("snGradlnk", this->runTime_.timeName(),this->mesh_),
        fvc::snGrad(lnk_)
    );
/*    volScalarField kDiffusion
    (
        "kDiffusion",
        fvc::laplacian(DkEff(), lnk_)
      + min(DkEff()*magSqr(fvc::grad(lnk_)), 1/dt)
//      + min(DkEff()*fvc::average(magSqr(fvc::snGrad(lnk_))), 1/dt)
    );
    kDiffusion.write();*/
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, lnk_)
      + fvm::div(alphaRhoPhi, lnk_)
      - fvm::laplacian(alpha*rho*DkEff(), lnk_)
      - min(DkEff()*magSqr(fvc::grad(lnk_)), 1/dt)
      //- min(DkEff()*fvc::average(magSqr(fvc::snGrad(lnk_))), 1/dt)
     ==
        alpha()*rho()*G/k1*exp(-lnk_())
      - 2.0/3.0*alpha()*rho()*divU
      - alpha()*rho()*eps1/k1*exp(lnepsilon_() - lnk_())
      //- this->nut_*N2_/(k1*sigmaTheta_)*exp(-lnk_)
      - Cmu_()*N2_()/sigmaTheta_*k1/eps1*exp(lnk_() - lnepsilon_())
      + fvModels.source(alpha, rho, lnk_)
    );

    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(lnk_);
    //bound(lnk_, log(this->kMin_/k1));
    k_ = k1*exp(lnk_);

    correctNut(gradU, S2, magS);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels

// Helper function
wordList swapWords
(
    const wordList& words,
    const word& oldWord,
    const word& newWord
)
{
    wordList newWords(words);
    
    forAll(newWords, i)
    {
        if (newWords[i] == oldWord) newWords[i] = newWord;
    }

    return newWords;
}

} // End namespace Foam

// ************************************************************************* //
