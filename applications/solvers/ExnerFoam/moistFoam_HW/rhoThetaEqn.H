// flux form rho theta equation

// Advect theta
dRhoThetadt = -fvc::div(U, theta);
theta = 
(
    rho.oldTime()*theta.oldTime() + dt*
    (
        (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
    )
)/rho;

// Calculate and apply source term
T = theta*Exner;
volScalarField heatSource("heatSource", theta*atmos.thetaSource(T, divu));

dRhoThetadt += rho*heatSource;
theta += dt*offCentre*heatSource;

// Update temperature
T = theta*Exner;

