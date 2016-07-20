// flux form rho theta equation

// Advect theta
dRhoThetadt = -fvc::div(U, theta);
theta = 
(
    theta.oldTime()*rho.oldTime() + dt*
    (
        (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
    )
)/rho;

// Calculate source term
volScalarField divu = fvc::div(U/rhof);
volScalarField cvml = Cv + Cvv*rv;
volScalarField cpml("cpml", Cp + Cpv*rv + rl*Cpl);
volScalarField Rm("Rm", R + Rv*rv);
T = TFromTheta(theta, Exner);
Lv = latentHeat(T, Lv0, T0, Cpl, Cpv);

volScalarField heatSource = theta*
(
  - (Rm/cvml - R*cpml/(Cp*cvml))*divu
  + (Cv*Lv/(cvml*Cp*T) - Rv/cvml*(1 - R/Cp*cpml/Rm))*condenseRate
);

// Apply source term
dRhoThetadt += rho*heatSource;
theta += dt*offCentre*heatSource;

// Update temperature
T = TFromTheta(theta, Exner);
