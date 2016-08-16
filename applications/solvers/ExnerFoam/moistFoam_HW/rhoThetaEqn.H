// flux form rho theta equation

// Advect thetaRho
dRhoThetadt = -fvc::div(un, rhoTheta);
rhoTheta = rhoTheta.oldTime() + dt*
(
    (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
);

// Calculate and apply source term
theta = rhoTheta/air.rho();
T = theta*Exner;
volScalarField heatSource("heatSource", rhoTheta*atmos.thetaSource(T, divu));

dRhoThetadt += heatSource;
rhoTheta += dt*offCentre*heatSource;

// Update temperature
theta = rhoTheta/air.rho();
T = theta*Exner;

//heatSource.write();
//divu.write();
//Exner.write();
//theta.write();
//atmos.write();

//water.updateCondensation(T);
//water.condensation().write();
//T.write();
//p.write();
//volScalarField ppGas("ppGas", water.gas().partialPressure(T));
//ppGas.write();
//volScalarField pSat("pSat", water.pSat(T));
//pSat.write();

//FatalErrorIn("") << exit(FatalError);

