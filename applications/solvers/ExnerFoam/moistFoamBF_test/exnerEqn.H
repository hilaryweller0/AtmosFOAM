exnerSource = -fvc::div(U,Exner)
            + Exner*divu
            - Exner*R*cpml/(Cp*cvml)*divu
            + R/Cp*(Lv/(cvml*theta) - Exner*Rv*cpml/(Rm*cvml))*condenseRate;

Exner = Exner.oldTime()
      + dt*(1-offCentre)*exnerSource.oldTime()
      + dt*offCentre*exnerSource;

