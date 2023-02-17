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

#include "stratifiedRealizableLogKE.H"
#include "fvOptions.H"
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

    return 1.0/(A0_ + As*Us*k1/eps1*exp(lnk_-lnepsilon_));
}


template<class BasicTurbulenceModel>
void stratifiedRealizableLogKE<BasicTurbulenceModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    Cmu_ = rCmu(gradU, S2, magS);
    this->nut_ = Cmu_*sqr(k1)/eps1*exp(2*lnk_ - lnepsilon_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void stratifiedRealizableLogKE<BasicTurbulenceModel>::correctNut()
{
    tmp<volTensorField> tgradU = fvc::grad(this->U_);

    surfaceTensorField gradUf = linearInterpolate(tgradU());
    gradUf += (fvc::snGrad(this->U_)*this->mesh_.magSf()
                 - (this->mesh_.Sf() & gradUf))
                *this->mesh_.Sf()/sqr(this->mesh_.magSf());
    surfaceScalarField S2f(2*magSqr(dev(symm(gradUf))));
    volScalarField S2 = fvc::average(S2f);

    //volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));
    correctNut(tgradU(), S2, magS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
stratifiedRealizableLogKE<BasicTurbulenceModel>::stratifiedRealizableLogKE
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
    lnk_
    (
        IOobject("lnk", this->runTime_.timeName(),this->mesh_),
        log(k_/k1),
        swapWords
        (
            k_.boundaryField().types(),
            "kqRWallFunction",
            "zeroGradient"
        )
    ),
    lnepsilon_
    (
        IOobject("lnepsilon", this->runTime_.timeName(),this->mesh_),
        log(epsilon_/eps1),
        swapWords
        (
            epsilon_.boundaryField().types(),
            "epsilonCmuWallFunction",
            "lnepsilonCmuWallFunction"
        )
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
    lnk_.writeOpt() = IOobject::AUTO_WRITE;
    lnepsilon_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool stratifiedRealizableLogKE<BasicTurbulenceModel>::read()
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
void stratifiedRealizableLogKE<BasicTurbulenceModel>::correct()
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
    fv::options& fvOptions(fv::options::New(mesh));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Intermediate fields
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    
    surfaceTensorField gradUf = linearInterpolate(tgradU());
    gradUf += (fvc::snGrad(U)*mesh.magSf() - (mesh.Sf() & gradUf))
                *mesh.Sf()/sqr(mesh.magSf());
    surfaceScalarField S2f(2*magSqr(dev(symm(gradUf))));
    volScalarField S2 = fvc::average(S2f);

    volScalarField magS(sqrt(S2));
    Cmu_ = rCmu(tgradU(), S2, magS);

    volScalarField eta(magS*k1/eps1*exp(lnk_ - lnepsilon_));
    volScalarField C1(max(eta/(scalar(5) + eta), scalar(0.43)));

    volScalarField G(this->GName(), nut*S2);

    // Brunt-Vaisala frequency squared
    N2_ = dbdz;
    bound(N2_, dimensionedScalar("", dbdz.dimensions(), scalar(0)));

    // Set ghat (gravity direction)
    dimensionedVector ghat("ghat", dimless, vector(0,0,1));

    // Richardson number
    volScalarField magSqrdUdz = magSqr(ghat & tgradU());
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
      + DepsilonEff()*magSqr(fvc::grad(lnepsilon_))
//      + DepsilonEff()*fvc::average(magSqr(fvc::snGrad(lnepsilon_)))
//      + DepsilonEff()*fvc::laplacian(lnepsilon_,lnepsilon_)
//      - DepsilonEff()*lnepsilon_*fvc::laplacian(lnepsilon_)
    );
    epsDiffusion.write();
*/
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, lnepsilon_)
      + fvm::div(alphaRhoPhi, lnepsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), lnepsilon_)
      - DepsilonEff()*magSqr(fvc::grad(lnepsilon_))
//      - alpha*rho*DepsilonEff()*fvc::laplacian(lnepsilon_, lnepsilon_)
//      + alpha*rho*DepsilonEff()*lnepsilon_*fvm::laplacian(lnepsilon_)
     ==
        C1*alpha*rho*magS
      - C2_*Cmu_*alpha*rho*k_/(this->nut_ + sqrt(this->nut_*this->nu()))
      - C3_*Cmu_*N2_/sigmaTheta_*k1/eps1*exp(lnk_-lnepsilon_)
      + C4_*min(sqrt(Ri_/C5_), scalar(1))*sqrt(N2_)
      + fvOptions(alpha, rho, lnepsilon_)
    );

    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(lnepsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(lnepsilon_);
    epsilon_ == eps1*exp(lnepsilon_);

    // Turbulent kinetic energy equation
    surfaceScalarField snGradlnk
    (
        IOobject("snGradlnk", this->runTime_.timeName(),this->mesh_),
        fvc::snGrad(lnk_)
    );
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, lnk_)
      + fvm::div(alphaRhoPhi, lnk_)
      - fvm::laplacian(alpha*rho*DkEff(), lnk_)
      - DkEff()*magSqr(fvc::grad(lnk_))
//      - alpha*rho*DkEff()*fvc::laplacian(lnk_, lnk_)
//      + alpha*rho*DkEff()*lnk_*fvm::laplacian(lnk_)
     ==
        alpha*rho*G/k1*exp(-lnk_)
      - 2.0/3.0*alpha*rho*divU
      - alpha*rho*eps1/k1*exp(lnepsilon_ - lnk_)
      //- this->nut_*N2_/(k1*sigmaTheta_)*exp(-lnk_)
      - Cmu_*N2_/sigmaTheta_*k1/eps1*exp(lnk_ - lnepsilon_)
      + fvOptions(alpha, rho, lnk_)
    );

    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(lnk_);
    k_ = k1*exp(lnk_);

    correctNut(tgradU(), S2, magS);
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
