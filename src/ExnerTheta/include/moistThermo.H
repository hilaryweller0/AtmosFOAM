/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Description
    Equation of State and pressure gradient coefficient for moist
    thermodynamics using Exner and thetaRho.

\*---------------------------------------------------------------------------*/

#ifndef moistThermo_H
#define moistThermo_H

#include "ExnerTheta.H"
#include "GeometricScalarField.H"
#include "surfaceFields.H"

#define TEMPLATE \
    template<class vector, template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

// The coefficient to convert the pressure gradient in terms of Exner to 
// grad(p)/rho
tmp<surfaceScalarField> gradPCoeff
(
    const volScalarField& cpml,
    const volScalarField& thetaRho,
    const volScalarField& rv,
    const scalar epsilon
)
{
    tmp<surfaceScalarField> tGradPCoeff
    (
        new surfaceScalarField
        (
            IOobject("gradPcoeff", thetaRho.time().timeName(), thetaRho.mesh()),
            fvc::interpolate(thetaRho*cpml/(1 + rv/epsilon))
        )
    );
    return tGradPCoeff;
}

// The coefficient to convert the pressure gradient in terms of Exner to 
// grad(p)/rho assuming dry Exner 
tmp<surfaceScalarField> gradPCoeff
(
    const volScalarField& theta,
    const volScalarField& rv,
    const volScalarField& rl,
    const dimensionedScalar& Cp,
    const scalar epsilon
)
{
    tmp<surfaceScalarField> tGradPCoeff
    (
        new surfaceScalarField
        (
            IOobject("gradPcoeff", theta.time().timeName(), theta.mesh()),
            Cp*fvc::interpolate(theta*(1+rv/epsilon)/(1+rv+rl))
        )
    );
    return tGradPCoeff;
}

// The pressure for the perfect gas law for a moist gas
tmp<volScalarField> moistGasP
(
    const volScalarField& rho,
    const volScalarField& Rm,
    const volScalarField& T,
    const volScalarField& rt
)
{
    tmp<volScalarField> tmoistGasP
    (
        new volScalarField
        (
            "p",
            rho*Rm*T/(1+rt)
        )
    );
    return tmoistGasP;
}

// The pressure for the perfect gas law for a moist gas
tmp<volScalarField> moistGasP
(
    const volScalarField& rho,
    const volScalarField& Rm,
    const volScalarField& T,
    const dimensionedScalar& rt0
)
{
    tmp<volScalarField> tmoistGasP
    (
        new volScalarField
        (
            "p",
            rho*Rm*T/(1+rt0)
        )
    );
    return tmoistGasP;
}

// The density in terms of p and T for the perfect gas law for a moist gas
tmp<volScalarField> moistGasPrimitiveRho
(
    const volScalarField& p,
    const volScalarField& Rm,
    const volScalarField& T,
    const volScalarField& rt
)
{
    tmp<volScalarField> tmoistGasPrimitiveRho
    (
        new volScalarField
        (
            "rho",
            p/(Rm*T)*(1+rt)
        )
    );
    return tmoistGasPrimitiveRho;
}

// The density in terms of p and T for the perfect gas law for a moist gas
tmp<volScalarField> moistGasPrimitiveRho
(
    const volScalarField& p,
    const volScalarField& Rm,
    const volScalarField& T,
    const dimensionedScalar& rt0
)
{
    tmp<volScalarField> tmoistGasPrimitiveRho
    (
        new volScalarField
        (
            "rho",
            p/(Rm*T)*(1+rt0)
        )
    );
    return tmoistGasPrimitiveRho;
}

// Exner in terms of rho and thetaRho for the perfect gas law for a moist gas
tmp<volScalarField> moistGasExner
(
    const volScalarField& rho,
    const volScalarField& thetaRho,
    const volScalarField& kappam,
    const dimensionedScalar& Rd,
    const dimensionedScalar& p0
)
{
    tmp<volScalarField> tmoistGasExner
    (
        new volScalarField
        (
            "Exner",
            pow(rho*Rd*thetaRho/p0, kappam/(1-kappam))
        )
    );
    return tmoistGasExner;
}

// Exner in terms of rho and theta for the perfect gas law for a moist gas
tmp<volScalarField> moistGasExner
(
    const volScalarField& rho,
    const volScalarField& theta,
    const volScalarField& rv,
    const dimensionedScalar& kappa,
    const dimensionedScalar& Rd,
    const dimensionedScalar& p0,
    const scalar epsilon
)
{
    tmp<volScalarField> tmoistGasExner
    (
        new volScalarField
        (
            "Exner",
            pow(rho*Rd*theta/p0*(1+rv/epsilon), kappa/(1-kappa))
        )
    );
    return tmoistGasExner;
}

// rho in terms of Exner and thetaRho for the perfect gas law for a moist gas
tmp<volScalarField> moistGasExnerRho
(
    const volScalarField& Exner,
    const volScalarField& thetaRho,
    const volScalarField& kappam,
    const dimensionedScalar& Rd,
    const dimensionedScalar& p0
)
{
    tmp<volScalarField> tmoistGasExnerRho
    (
        new volScalarField
        (
            "rho",
            p0/(Rd*thetaRho)*pow(Exner, (1-kappam)/kappam)
        )
    );
    return tmoistGasExnerRho;
}

// rho in terms of Exner and thetaRho for the perfect gas law for a moist gas
tmp<volScalarField> moistGasExnerRho
(
    const volScalarField& Exner,
    const volScalarField& theta,
    const volScalarField& rv,
    const dimensionedScalar& kappa,
    const dimensionedScalar& Rd,
    const dimensionedScalar& p0,
    const scalar epsilon
)
{
    tmp<volScalarField> tmoistGasExnerRho
    (
        new volScalarField
        (
            "rho",
            p0/(Rd*theta)*pow(Exner, (1-kappa)/kappa)/(1+rv/epsilon)
        )
    );
    return tmoistGasExnerRho;
}

// Primitive pressure, p, from Exner
tmp<volScalarField> pFromExner
(
    const volScalarField& Exner,
    const volScalarField& kappam,
    const dimensionedScalar& p0
)
{
    tmp<volScalarField> tpFromEnxer
    (
        new volScalarField
        (
            "p",
            p0*pow(Exner, 1/kappam)
        )
    );
    return tpFromEnxer;
}

// Primitive pressure, p, from Exner
tmp<volScalarField> pFromExner
(
    const volScalarField& Exner,
    const dimensionedScalar& kappa,
    const dimensionedScalar& p0
)
{
    tmp<volScalarField> tpFromEnxer
    (
        new volScalarField
        (
            "p",
            p0*pow(Exner, 1/kappa)
        )
    );
    return tpFromEnxer;
}

// Primitive temperature, T, from thetaRho
tmp<volScalarField> TFromThetaRho
(
    const volScalarField& thetaRho,
    const volScalarField& Exner,
    const volScalarField& rv,
    const volScalarField& rl,
    const scalar epsilon
)
{
    tmp<volScalarField> tTFromThetaRho
    (
        new volScalarField
        (
            "T",
            thetaRho*Exner*(1+rv+rl)/(1+rv/epsilon)
        )
    );
    return tTFromThetaRho;
}

// Primitive temperature, T, from thetaRho
tmp<volScalarField> TFromThetaRho
(
    const volScalarField& thetaRho,
    const volScalarField& Exner,
    const volScalarField& rv,
    const dimensionedScalar& rt0,
    const scalar epsilon
)
{
    tmp<volScalarField> tTFromThetaRho
    (
        new volScalarField
        (
            "T",
            thetaRho*Exner*(1+rt0)/(1+rv/epsilon)
        )
    );
    return tTFromThetaRho;
}

// Primitive temperature, T, from theta
tmp<volScalarField> TFromTheta
(
    const volScalarField& theta,
    const volScalarField& Exner
)
{
    tmp<volScalarField> tTFromTheta
    (
        new volScalarField
        (
            "T",
            theta*Exner
        )
    );
    return tTFromTheta;
}

// thetaRho from temperature, T,
tmp<volScalarField> thetaRhoFromT
(
    const volScalarField& T,
    const volScalarField& Exner,
    const volScalarField& rv,
    const dimensionedScalar& rt0,
    const scalar epsilon
)
{
    tmp<volScalarField> tthetaRhoFromT
    (
        new volScalarField
        (
            "thetaRho",
            T/Exner*(1 + rv/epsilon)/(1+rt0)
        )
    );
    return tthetaRhoFromT;
}

// theta from temperature, T,
tmp<volScalarField> thetaRhoFromT
(
    const volScalarField& T,
    const volScalarField& Exner
)
{
    tmp<volScalarField> tthetaFromT
    (
        new volScalarField
        (
            "theta",
            T/Exner
        )
    );
    return tthetaFromT;
}

// Latent heat of vaporisation
tmp<volScalarField> latentHeat
(
    const volScalarField& T,
    const dimensionedScalar& Lv0,
    const dimensionedScalar& T0,
    const dimensionedScalar& Cpl,
    const dimensionedScalar& Cpv
)
{
    tmp<volScalarField> tlatentHeat
    (
        new volScalarField
        (
            "Lv",
            Lv0 - (T-T0)*(Cpl - Cpv)
        )
    );
    return tlatentHeat;
}

// Saturation vapour pressure from the Claussiua Clapyron relation
tmp<volScalarField> pSat
(
    const volScalarField& T,
    const dimensionedScalar& pSat0,
    const dimensionedScalar& Lv0,
    const dimensionedScalar& Rv,
    const dimensionedScalar& T0
)
{
    tmp<volScalarField> tes
    (
        new volScalarField
        (
            "es",
            pSat0*Foam::exp(-Lv0/Rv*(1/T - 1/T0))
        )
    );
    return tes;
}

// Equivalent potential temperature from temperature and pressure
volScalarField thetaeFromPrimitive
(
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& Lv,
    const volScalarField& rv,
    const volScalarField& rl,
    const dimensionedScalar& p0,
    const dimensionedScalar& epsilon,
    const dimensionedScalar& Cpd,
    const dimensionedScalar& Cpl,
    const dimensionedScalar& Rd
)
{
    volScalarField Cp = (rl+rv)*Cpl + Cpd;

    return T*pow(p/p0*epsilon/(rv+epsilon), -Rd/Cp)
           *exp(Lv*rv/(Cp*T));
}


// Equivalent potential temperature from temperature and pressure
tmp<volScalarField> thetaeFromPrimitive
(
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& Lv,
    const volScalarField& rv,
    const dimensionedScalar& rt0,
    const dimensionedScalar& p0,
    const dimensionedScalar& epsilon,
    const dimensionedScalar& Cpd,
    const dimensionedScalar& Cpl,
    const dimensionedScalar& Rd
)
{
    dimensionedScalar Cp = Cpd + Cpl*rt0;

    tmp<volScalarField> tthetae
    (
        new volScalarField
        (
            "thetae",
            T*pow(p/p0/(1+rv/epsilon), -Rd/Cp)
           *exp(Lv*rv/(Cp*T))
        )
    );
    return tthetae;
}


// Iterations to find temperature, T, from Equivalent potential temperature
void TfromThetae
(
    volScalarField& T,
    const dimensionedScalar& thetae0,
    const volScalarField& p,
    const volScalarField& Lv,
    const volScalarField& rv,
    const dimensionedScalar& rt0,
    const dimensionedScalar& p0,
    const dimensionedScalar& epsilon,
    const dimensionedScalar& Cpd,
    const dimensionedScalar& Cpl,
    const dimensionedScalar& Rd
)
{
    // Setup variables for a Newton method
    dimensionedScalar Cp = Cpd + Cpl*rt0;
    volScalarField a = pow(p/p0/(1+rv/epsilon), -Rd/Cp);
    volScalarField b = Lv*rv/Cp;

    // The residual to miniminse every Newton step
    volScalarField resid = thetae0 - T*a*exp(b/T);
    scalar RMSresid = sqrt(sum(sqr(resid.internalField())));
    
    // convergence criterial
    scalar conv = 1e-6;
    
    // Newton steps
    for(label i = 0; i < 10 && RMSresid > conv; i++)
    {
        T -= resid/(a*exp(b/T)*(b/T - 1));
        resid = thetae0 - T*a*exp(b/T);
        RMSresid = sqrt(sum(sqr(resid.internalField())));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#endif
// ************************************************************************* //
