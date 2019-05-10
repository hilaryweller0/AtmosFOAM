// advective form theta equation
dThetaDt = -fvc::div(volFlux, theta, "theta") + theta*fvc::div(volFlux);

// Viscous term
if (alpha.value() > VSMALL)
{
    dThetaDt += (1/rho)*fvc::laplacian(turbulence->alphaEff(), theta);
}

// Update theta according to the rate of change
theta = theta.oldTime() + dt*
(
    (1-offCentre)*dThetaDt.oldTime()
  + offCentre*dThetaDt
);

theta.correctBoundaryConditions();
thetaf = fvc::interpolate(theta);

thermo.T() = theta*Exner;
