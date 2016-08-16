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

