// flux form rho theta equation

// Variables needed to calculate the moist source term
cpml = Cp + Cpv*rv;
Rm = R + Rv*rv;
Lv = Lv0 - (Cpl - Cpv)*(T - T0);

dRhoThetadt = thetaRho/qd*
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

T = thetaRho*Exner/(1 + qv*(Rv/R-1) - ql);

