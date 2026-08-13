/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "stableBLkEpsilon.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> stableBLkEpsilon<BasicMomentumTransportModel>::boundEpsilon()
{
    tmp<volScalarField> tCmuk2(this->Cmu_*sqr(this->k_));
    this->epsilon_ = max(this->epsilon_, this->Cmu_*sqr(this->kMin_)/this->nu());
    return tCmuk2;
}


template<class BasicMomentumTransportModel>
void stableBLkEpsilon<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = boundEpsilon()/this->epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
stableBLkEpsilon<BasicMomentumTransportModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        return fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return kEpsilon<BasicMomentumTransportModel>::kSource();
    }
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
stableBLkEpsilon<BasicMomentumTransportModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        volScalarField Gpos = max
        (
            Gcoef(),
            dimensionedScalar(dimDensity/dimTime, scalar(0))
        );

        return fvm::SuSp(this->C1_*Gpos, this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicMomentumTransportModel>::epsilonSource();
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
stableBLkEpsilon<BasicMomentumTransportModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    return (Cg_*this->Cmu_)*this->alpha_*this->rho_*this->k_
           *(g & fvc::grad(theta_))/(theta_*this->epsilon_);
}


template<class BasicMomentumTransportModel>
void stableBLkEpsilon<BasicMomentumTransportModel>::updateSurfaceFields()
{
    // Look up fields for the surface
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
    const volVectorField& U = 
        this->mesh_.objectRegistry::template
        lookupObject<volVectorField>("U");

    volScalarField::Boundary& uStarBf = uStar_.boundaryFieldRef();
    volScalarField::Boundary& thetaStarBf = thetaStar_.boundaryFieldRef();
    volScalarField::Boundary& LmoBf = Lmo_.boundaryFieldRef();
    volScalarField::Boundary& nutG = nutGround_.boundaryFieldRef();
    volScalarField::Boundary& alphatG = alphatGround_.boundaryFieldRef();
    volScalarField::Boundary& epsilonGcorr
         = epsilonGroundCorr_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U.boundaryField();
    const volScalarField::Boundary& thetaBf = theta_.boundaryField();
    const volScalarField::Boundary& z0Bf = z0_.boundaryField();

    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& y = this->yb()[patchi];
            //uStarBf[patchi] = sqrt(mag(devTauBf[patchi]));
            forAll(thetaStarBf[patchi], facei)
            {
                const label celli = patch.faceCells()[facei];
                
                uStarBf[patchi][facei] = kappa_.value()*
                (
                    mag(U[celli]) - mag(UBf[patchi][facei])
                )/
                (
                    log(y[facei]/max(z0Bf[patchi][facei],SMALL))
                  + betam_.value()*y[facei]/LmoBf[patchi][facei]
                );
                
                thetaStarBf[patchi][facei] = max
                (
                    kappa_.value()*(theta_[celli] - thetaBf[patchi][facei])/
                    (
                        log(y[facei]/max(z0Bf[patchi][facei],SMALL))
                      + betah_.value()*y[facei]/LmoBf[patchi][facei]
                    ),
                    SMALL
                );
            }
            LmoBf[patchi] = Tref_.value()*sqr(uStarBf[patchi])
                         /(mag(g.value())*thetaStarBf[patchi]*kappa_.value());
            
            nutG[patchi] = uStarBf[patchi]*y*kappa_.value()/
            (
                log(y/max(z0Bf[patchi],SMALL))
              + betam_.value()*y/LmoBf[patchi]
            );
            alphatG[patchi] = thetaStarBf[patchi]*y*kappa_.value()/
            (
                log(y/max(z0Bf[patchi],SMALL))
              + betah_.value()*y/LmoBf[patchi]
            );
            epsilonGcorr[patchi] = (1 + (betam_.value()-1)*y/LmoBf[patchi]);
        }
    }
}


template<class BasicMomentumTransportModel>
volScalarField& stableBLkEpsilon<BasicMomentumTransportModel>::calcTheta()
{
    const volScalarField& T =
        this->mesh_.objectRegistry::template lookupObject<volScalarField>("T");
    const volScalarField& p =
        this->mesh_.objectRegistry::template lookupObject<volScalarField>("p");

    theta_ = T*pow(pRef_/p, RbyCp_);

    return theta_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
stableBLkEpsilon<BasicMomentumTransportModel>::stableBLkEpsilon
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
    kEpsilon<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),
    Cg_("Cg", this->typeDict(type), 1.0),
    kappa_("kappa", this->typeDict(type), 0.41),
    betam_("betam", this->typeDict(type), 4.8),
    betah_("betah", this->typeDict(type), 7.8),
    Tref_("Tref",   dimTemperature, this->typeDict(type)),
    z0_
    (
        IOobject("z0", "constant", this->mesh_, IOobject::MUST_READ),
        this->mesh_
    ),
    Lmo_
    (
        IOobject("Lmo", this->runTime_.name(), this->mesh_, 
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        this->mesh_,
        dimensionedScalar(dimLength, GREAT)
    ),
    uStar_
    (
        IOobject("uStar", this->runTime_.name(), this->mesh_, 
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        this->mesh_,
        dimensionedScalar(dimVelocity, scalar(0))
    ),
    thetaStar_
    (
        IOobject("thetaStar", this->runTime_.name(), this->mesh_,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        this->mesh_,
        dimensionedScalar(dimTemperature, scalar(0))
    ),
    nutGround_
    (
        IOobject("nutGround", this->runTime_.name(), this->mesh_),
        this->mesh_,
        dimensionedScalar(dimensionSet(0,2,-1,0,0), scalar(0))
    ),
    alphatGround_
    (
        IOobject("alphatGround", this->runTime_.name(), this->mesh_),
        this->mesh_,
        dimensionedScalar(dimensionSet(0,2,-1,0,0), scalar(0))
    ),
    epsilonGroundCorr_
    (
        IOobject("epsilonGroundCorr", this->runTime_.name(), this->mesh_),
        this->mesh_,
        dimensionedScalar(dimless, scalar(0))
    ),
    theta_
    (
        IOobject("theta", this->runTime_.name(), this->mesh_,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        this->mesh_,
        Tref_
    ),
    pRef_
    (
        this->mesh_.objectRegistry::template
            lookupObject<uniformDimensionedScalarField>("pRef")
    ),
    RbyCp_(0)
{
    IOdictionary thermoDict
    (
        this->mesh_.objectRegistry::template 
            lookupObject<IOdictionary>("physicalProperties")
    );
    scalar W = readScalar
    (
        thermoDict.subDict("mixture").subDict("specie").lookup("molWeight")
    );
    scalar Cp = readScalar
    (
        thermoDict.subDict("mixture").subDict("thermodynamics").lookup("Cp")
    );
    RbyCp_ = constant::thermodynamic::RR/(W*Cp);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool stableBLkEpsilon<BasicMomentumTransportModel>::read()
{
    if (kEpsilon<BasicMomentumTransportModel>::read())
    {
        Cg_.readIfPresent(this->typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void stableBLkEpsilon<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    calcTheta();
    updateSurfaceFields();
    kEpsilon<BasicMomentumTransportModel>::correct();
    updateSurfaceFields();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
