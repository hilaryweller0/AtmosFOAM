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
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicTurbulenceModel>::kSource() const
{
    return fvm::SuSp(-this->nut_/this->k_*N2_/sigmaTheta_, this->k_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> stratifiedRealizableKE<BasicTurbulenceModel>::epsilonSource() const
{
    return fvm::SuSp
    (
        C3_/this->k_*this->nut_*N2_/sigmaTheta_
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
    realizableKE<BasicTurbulenceModel>
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
    sigmaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaTheta",
            this->coeffDict_,
            0.74
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
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool stratifiedRealizableKE<BasicTurbulenceModel>::read()
{
    if (realizableKE<BasicTurbulenceModel>::read())
    {
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
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
    const volVectorField& U = this->U_;

    // Fields for calculating Brunt-Vaisala and Richardson no.
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
    const dimensionedVector ghat = g/mag(g);

    // Define Brunt-Vaisala frequency and Richardson number depending on the 
    // data avaialble. Use one of
    // gradThetaByTheta, theta, dbdz or (p and T)
    if
    (
        this->mesh_.objectRegistry::template
             foundObject<volVectorField>("gradThetaByTheta")
    )
    {
        const volVectorField& gradThetaByTheta
             = this->mesh_.objectRegistry::template
                 lookupObject<volVectorField>("gradThetaByTheta");
    
        // Brunt-Vaisala frequency squared
        N2_ = g & gradThetaByTheta;

        // Richardson number
        Ri_ = N2_/magSqr(fvc::grad(U) & ghat);
    }
    else if
    (
        this->mesh_.objectRegistry::template
             foundObject<volScalarField>("theta")
    )
    {
        const volScalarField& theta = this->mesh_.objectRegistry::template
                 lookupObject<volScalarField>("theta");
        const volVectorField gradThetaByTheta = fvc::grad(theta)/theta;
    
        // Brunt-Vaisala frequency squared
        N2_ = g & gradThetaByTheta;

        // Richardson number
        Ri_ = N2_/magSqr(fvc::grad(U) & ghat);
    }
    else if
    (
        this->mesh_.objectRegistry::template foundObject<volScalarField>("dbdz")
    )
    {
        const volScalarField& dbdz = this->mesh_.objectRegistry::template
                 lookupObject<volScalarField>("dbdz");
    
        // Brunt-Vaisala frequency squared
        N2_ = dbdz;

        // Richardson number
        Ri_ = N2_/magSqr(fvc::grad(U) & ghat);
    }
    else // Need to calculate theta from T and p from the transport model
    {
        const volScalarField& T = this->transport_.T();
        const volScalarField& p = this->transport_.p();

        // Calculate potential temperature
        const volScalarField kappa
        (
            "kappa",
             1-1/this->transport_.gamma()
        );
        const dimensionedScalar pRef("pRef", dimPressure, 1e5);
        volScalarField theta = T*pow(pRef/p, kappa);
        // Smooth theta
        //surfaceScalarField thetaf = linearInterpolate(theta);
        //theta = fvc::average(thetaf);

        // Brunt-Vaisala frequency squared
        N2_ = (g & (fvc::grad(theta)))/theta;

        // Richardson number
        Ri_ = N2_/magSqr(fvc::grad(U) & ghat);
    }
    
    realizableKE<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
