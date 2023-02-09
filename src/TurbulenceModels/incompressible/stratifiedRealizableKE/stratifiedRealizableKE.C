/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
#include "fvOptions.H"
#include "bound.H"
#include "uniformDimensionedFields.H"

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

    volScalarField W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
        )
    );

    tS.clear();

    volScalarField phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    volScalarField As(sqrt(6.0)*cos(phis));
    volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    /*dimensionedScalar epsSMALL("epsSMALL", epsilon_.dimensions(), 1e8*SMALL);
    Info << "k = " << k_[999] << " A0 = " << A0_.value() << " As = " << As[999]
         << " Us = " << Us[999] << " epsilon = " << epsilon_[999] 
         << " epsSMALL = " << epsSMALL.value()
         << " kMin = " << this->kMin_.value()
         << " epsilonMin = " << this->epsilonMin_.value()
         << " VSMALL = " << VSMALL   << endl;
    volScalarField Cmu1 = 1.0/(A0_ + As*Us*k_/epsilon_);
    volScalarField Cmu2 = (epsilon_+0.09*A0_*epsSMALL)
                        /(A0_*(epsilon_+epsSMALL) + As*Us*k_);
    Info << "Cmu1 = " << Cmu1[999] << " Cmu2 = " << Cmu2[999] << endl;*/
    //return (epsilon_+0.09*A0_*epsSMALL)/(A0_*(epsilon_+epsSMALL) + As*Us*k_);
    return 1.0/(A0_ + As*Us*k_/epsilon_);
}


template<class BasicTurbulenceModel>
void stratifiedRealizableKE<BasicTurbulenceModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    Cmu_ = rCmu(gradU, S2, magS);
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void stratifiedRealizableKE<BasicTurbulenceModel>::correctNut()
{
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));
    correctNut(tgradU(), S2, magS);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicTurbulenceModel>::kSource() const
{
    return fvm::SuSp
    (
        -Cmu_*N2_*this->k_/(this->epsilon_*sigmaTheta_),
        this->k_
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicTurbulenceModel>::epsilonSource() const
{
    return fvm::SuSp
    (
        -C3_*Cmu_*this->k_*N2_/(this->epsilon_*sigmaTheta_)
      + C4_*min(sqrt(Ri_/C5_), scalar(1))*sqrt(N2_),
        this->epsilon_
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
stratifiedRealizableKE<BasicTurbulenceModel>::stratifiedRealizableKE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
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
    writeCmu_
    (
        Switch::lookupOrAddToDict("writeCmu", this->coeffDict_, false)
    ),
    k_
    (
        IOobject
        (
            "k",
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
            "epsilon",
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
        dimensionedScalar("Cmu", dimless, scalar(0))
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    if (writeCmu_)
    {
        Cmu_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool stratifiedRealizableKE<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
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


template<class BasicTurbulenceModel>
void stratifiedRealizableKE<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Calculate turbulence parameters not based on stratification
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));

    volScalarField eta(magS*k_/epsilon_);
    volScalarField C1(max(eta/(scalar(5) + eta), scalar(0.43)));

    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));

    // Set ghat (gravity direction)
    dimensionedVector ghat("ghat", dimless, vector(0,0,1));

    // Define Brunt-Vaisala frequency
    const volScalarField& dbdz = this->mesh_.objectRegistry::template
                 lookupObject<volScalarField>("dbdz");
    // Brunt-Vaisala frequency squared
    N2_ = dbdz;

    // Ensure that N2 is not negative
    N2_ = max(N2_, dimensionedScalar("", N2_.dimensions(), scalar(0)));

    // Richardson number
    Ri_ = N2_/max(magSqr(ghat & tgradU()), this->epsilonMin()/this->nu());
    
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1*alpha*rho*magS*epsilon_
      - fvm::Sp
        (
            C2_*alpha*rho*epsilon_/(k_ + sqrt(this->nu()*epsilon_)),
            epsilon_
        )
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilon_/k_, k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(tgradU(), S2, magS);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
