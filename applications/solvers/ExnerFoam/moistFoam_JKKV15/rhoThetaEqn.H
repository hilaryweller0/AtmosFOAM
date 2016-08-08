// flux form rho theta equation

// Advect thetaRho
dRhoThetadt = -fvc::div(U, thetaRho);
thetaRho = 
(
    thetaRho.oldTime()*rho.oldTime() + dt*
    (
        (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
    )
)/rho;

// Calculate and apply source term
cpml = Cp + rv*Cpv + rl*Cpl;
Rm = R + Rv*rv;
T = TFromThetaRho(thetaRho, Exner, rv, rl, epsilon);
Lv = latentHeat(T, Lv0, T0, Cpl, Cpv);

volScalarField heatSource = -thetaRho/qd*
(
    Rv/Rm
  - Foam::log(Exner)*(Rv/Rm + (Cpl - Cpv)/cpml)
  - Lv/(cpml*T)
)*condenseRate;

dRhoThetadt += rho*heatSource;
thetaRho += dt*offCentre*heatSource;

// Update temperature
T = TFromThetaRho(thetaRho, Exner, rv, rl, epsilon);
