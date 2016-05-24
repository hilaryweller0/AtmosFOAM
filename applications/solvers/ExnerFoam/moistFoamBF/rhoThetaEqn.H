// flux form rho theta equation

// Variables needed to calculate the moist source term
volScalarField divu = fvc::div(U/rhof);
volScalarField cvml = Cv + Cvv*rv;
volScalarField cpml = Cp + Cpv*rv;
volScalarField Rm = R + Rv*rv;

rhoThetaSource = rho*theta*
(
  - (Rm/cvml - R*cpml/(Cp*cvml))*divu
  + (Cv*Lv/(cvml*Cp*T) - Rv/cvml*(1 - R/Cp*cpml/Rm))*condenseRate
);

theta = 
(
    theta.oldTime()*rho.oldTime() - dt*
    (
        (1-offCentre)*divUtheta.oldTime()
      + offCentre*fvc::div(U,theta)
      - (1-offCentre)*rhoThetaSource.oldTime()
      - offCentre*rhoThetaSource
    )
)/rho;

T = theta*Exner;

