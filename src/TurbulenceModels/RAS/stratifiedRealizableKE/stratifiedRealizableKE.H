/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::RASModels::stratifiedRealizableKE

Description
    k and epsilon model suitable for neutral and stratified flow from

    Reference:
    \verbatim
        Xiping Zeng, Yansen Wang, and Benjamin T. MacCall. (2020)
        A k − ε turbulence model for the stable atmosphere
        J. of Atmos Sci., 77(1):167–184, 2020
    \endverbatim

    The model coefficients from ZWM20 are
    \verbatim
        stratifiedRealizableKECoeffs
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
    stratifiedRealizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef stratifiedRealizableKE_H
#define stratifiedRealizableKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class stratifiedRealizableKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class stratifiedRealizableKE
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
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
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("stratifiedRealizableKE");


    // Constructors

        //- Construct from components
        stratifiedRealizableKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );


    //- Destructor
    virtual ~stratifiedRealizableKE()
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

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(0.09*k_),
                epsilon_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stratifiedRealizableKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
