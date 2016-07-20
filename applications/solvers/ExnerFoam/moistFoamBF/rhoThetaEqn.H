// flux form rho theta equation

// Variables needed to calculate the moist source term
volScalarField divu = fvc::div(U/rhof);
volScalarField cvml = Cv + Cvv*rv;
volScalarField cpml("cpml", Cp + Cpv*rv + rl*Cpl);
volScalarField Rm("Rm", R + Rv*rv);
T = TFromTheta(theta, Exner);
Lv = latentHeat(T, Lv0, T0, Cpl, Cpv);

// Could put splitting into theta eqn

dRhoThetadt = theta*
(
  - (Rm/cvml - R*cpml/(Cp*cvml))*divu*rho
  + (Cv*Lv/(cvml*Cp*T) - Rv/cvml*(1 - R/Cp*cpml/Rm))*condenseRate
)
- fvc::div(U,theta);

theta = 
(
    theta.oldTime()*rho.oldTime() + dt*
    (
        (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
    )
)/rho;

T = TFromTheta(theta, Exner);
