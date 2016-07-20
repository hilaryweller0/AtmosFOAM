// flux form rho theta equation

// Variables needed to calculate the moist source term
cpml = Cp + rv*Cpv + rl*Cpl;
Rm = R + Rv*rv;
T = TFromThetaRho(thetaRho, Exner, rv, rl, epsilon);
Lv = latentHeat(T, Lv0, T0, Cpl, Cpv);

dRhoThetadt = -thetaRho/qd*
(
    Rv/Rm
  - Foam::log(Exner)*(Rv/Rm + (Cpl - Cpv)/cpml)
  - Lv/(cpml*T)
)*condenseRate
- fvc::div(U, thetaRho);

thetaRho = 
(
    thetaRho.oldTime()*rho.oldTime() + dt*
    (
        (1-offCentre)*dRhoThetadt.oldTime() + offCentre*dRhoThetadt
    )
)/rho;

T = TFromThetaRho(thetaRho, Exner, rv, rl, epsilon);
