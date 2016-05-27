// flux form rho theta equation

// Variables needed to calculate the moist source term
volScalarField cvml = Cv + Cvv*rv;
volScalarField cpml = Cp + Cpv*rv;
volScalarField Rm = R + Rv*rv;

thetaSource = -fvc::div(U,theta)
            + theta*divu
            - theta*(Rm/cvml - R*cpml/(Cp*cvml))*divu
            + (
                  Cv*Lv/(cvml*Cp*Exner)
                - theta*Rv/cvml*(1 - R/Cp*cpml/Rm)
              )*condenseRate;

theta = theta.oldTime()
      + dt*(1-offCentre)*thetaSource.oldTime()
      + dt*offCentre*thetaSource;

