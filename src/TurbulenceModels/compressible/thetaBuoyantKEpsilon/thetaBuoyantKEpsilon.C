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

#include "thetaBuoyantKEpsilon.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
thetaBuoyantKEpsilon<BasicTurbulenceModel>::thetaBuoyantKEpsilon
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
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool thetaBuoyantKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
thetaBuoyantKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    const volScalarField& T = this->transport_.T();
    const volScalarField& p = this->transport_.p();
    const volScalarField& rho = this->rho_;
    
    const volScalarField kappa
    (
        "kappa",
         1-1/this->transport_.gamma()
    );
    const dimensionedScalar pRef("pRef", dimPressure, 1e5);
    volScalarField theta = T*pow(pRef/p, kappa);
    surfaceScalarField thetaf = linearInterpolate(theta);
    theta = fvc::average(thetaf);

    tmp<volScalarField> Gcoefp
    (
        new volScalarField
        (
            "Gcoef", 
            -(Cg_*this->Cmu_)*this->alpha_*rho*this->k_/kappa*
            (
                g &
                (
                    fvc::grad(theta)/theta
                )
            )
           /(this->epsilon_ + this->epsilonMin_)
        )
    );

    return Gcoefp;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
thetaBuoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        return -fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::kSource();
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
thetaBuoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        vector gHat(g.value()/mag(g.value()));

        volScalarField w(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*w)
          + dimensionedScalar(dimVelocity, small)
        );

        return -fvm::SuSp(this->C1_*tanh(mag(w)/u)*Gcoef(), this->epsilon_);

//        return -fvm::SuSp(this->C1_*Gcoef(), this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::epsilonSource();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
