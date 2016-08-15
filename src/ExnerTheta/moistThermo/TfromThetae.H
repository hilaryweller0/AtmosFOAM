// Newton iterations to find temperature, T, from equivalent potential
// temperature, thetae, given a uniform total mixing ratio, rt

#include "atmosphere.H"

void TfromThetae
(
    volScalarField& T,
    atmosphere& atmos,
    const dimensionedScalar& thetae0,
    const scalar rt
)
{
    perfectGasPhase& air = atmos[0].gas();
    fluidSpecie& water = atmos[1];

    // Setup variables for a Newton method
    dimensionedScalar Cp = air.Cp() + water.liquid().Cp()*rt;
    volScalarField a = pow(air.partialPressure(T)/air.p0(), -air.R()/Cp);
    volScalarField b = water.gas().rho()*water.latentHeat(T)/(Cp*air.rho());

    // The residual to miniminse every Newton step
    volScalarField resid = thetae0 - T*a*exp(b/T);
    scalar RMSresid = Foam::sqrt(sum(sqr(resid.internalField())).value());
    
    // convergence criterial
    scalar conv = 1e-6;
    
    // Newton steps
    for(label i = 0; i < 10 && RMSresid > conv; i++)
    {
        T -= resid*T/(a*(b-T)*exp(b/T));
        resid = thetae0 - T*a*exp(b/T);
        RMSresid = Foam::sqrt(sum(sqr(resid.internalField())).value());
    }
}
