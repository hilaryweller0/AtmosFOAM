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

Class
    Foam::RASModels::stratifiedRealizableLogKE

Description
    k and epsilon model suitable for neutral and stratified flow from

    Reference:
    \verbatim
        Xiping Zeng, Yansen Wang, and Benjamin T. MacCall. (2020)
        A k − ε turbulence model for the stable atmosphere
        J. of Atmos Sci., 77(1):167–184, 2020
    \endverbatim
    Solves for log k and log epsilon

    The model coefficients from ZWM20 are
    \verbatim
        stratifiedRealizableLogKECoeffs
        {
            A0              4.0;
            C2              1.9;
            C3              1.46;
            C4              0.44;
            C5              0.08;
            sigmak          1;
            sigmaEps        1.2;
            sigmaTheta      0.74;
        }
    \endverbatim
    Other variables are
            C1 = max(0.43, eta/(eta+5);
            Cmu = 1/(A0 + k/eps sqrt(6(S:S + Omega:Omega)) cos(phi))
            eta = meanS k/eps;
            Omega is the mean rotation rate tensor ...
            S = 0.5(gradU + gradU^T)
            meanS = sqrt(2 S:S)

See also
    Foam::RASModels::realizableKE

SourceFiles
    stratifiedRealizableLogKE.C

\*---------------------------------------------------------------------------*/

#ifndef stratifiedRealizableLogKE_H
#define stratifiedRealizableLogKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

static dimensionedScalar k1
(
    "k1",
    dimensionSet(0,2,-2,0,0),
    scalar(1)
);

static dimensionedScalar eps1
(
    "eps1",
    dimensionSet(0,2,-3,0,0),
    scalar(1)
);

/*---------------------------------------------------------------------------*\
                           Class stratifiedRealizableLogKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class stratifiedRealizableLogKE
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaTheta_;

        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField lnk_;
            volScalarField lnepsilon_;

            // Square of the Brunt-Vaisala frequency
            volScalarField N2_;
            // Bulk Richardson number
            volScalarField Ri_;
            // Cmu
            volScalarField Cmu_;

   // Protected Member Functions

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("stratifiedRealizableLogKE");


    // Constructors

        //- Construct from components
        stratifiedRealizableLogKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~stratifiedRealizableLogKE()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the effective diffusivity for theta
        tmp<volScalarField> DthetaEff() const
        {
            return volScalarField::New
            (
                "DthetaEff",
                ((this->nut_ + this->nu())/sigmaTheta_)
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return Cmu, Ri and N2
        tmp<volScalarField> Cmu() const { return Cmu_; }
        tmp<volScalarField> Ri() const { return Ri_; }
        tmp<volScalarField> N2() const { return N2_; }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stratifiedRealizableLogKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
