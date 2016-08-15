// flux form rho theta equation

// Advect thetaRho
dRhoThetadt = -fvc::div(un, rhoTheta);
rhoTheta = rhoTheta.oldTime() + dt*
(
    (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
);

// Calculate and apply source term
theta = rhoTheta/atmos.rho();
T = theta*Exner;
volScalarField heatSource = rhoTheta*atmos.thetaSource(T, divu);

dRhoThetadt += heatSource;
rhoTheta += dt*offCentre*heatSource;

// Update temperature
theta = rhoTheta/atmos.rho();
T = theta*Exner;

