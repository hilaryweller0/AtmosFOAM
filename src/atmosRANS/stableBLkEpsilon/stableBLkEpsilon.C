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
        volScalarField Gneg = min
        (
            Gcoef(),
            dimensionedScalar(dimDensity/dimTime, scalar(0))
        );

        return fvm::SuSp(this->C1_*Gneg, this->epsilon_);
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

//    const uniformDimensionedScalarField& pRef =
//        this->mesh_.objectRegistry::template
//        lookupObject<uniformDimensionedScalarField>("pRef");
    
    const volScalarField& T = 
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>("T");

    return
        (Cg_*this->Cmu_)*this->alpha_*this->rho_*this->k_*(g & fvc::grad(T))/T
       /this->epsilon_;
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

    Cg_("Cg", this->typeDict(type), 1.0)
{}


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



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
