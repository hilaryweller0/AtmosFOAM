// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", Pcc*pRef*Foam::exp(-Lv0/Rv*(1/T - 1/T0)));

// Calculate the saturation mixing ratio
volScalarField rvs = epsilon*es/(p-es);
Lv = Lv0 - (Cpl - Cpv)*(T - T0);

// Update rv and rl without splitting
rv =
(
    rv.oldTime()*rho.oldTime() - dt*
    (
        (1-offCentre)*divUrv.oldTime()
      + offCentre*fvc::div(U,rv)
      + (1-offCentre)*condenseRate.oldTime()*rho.oldTime()
      + offCentre*condenseRate*rho
    )
)/rho;

rl =
(
    rl.oldTime()*rho.oldTime() - dt*
    (
        (1-offCentre)*divUrl.oldTime()
      + offCentre*fvc::div(U,rl)
      - (1-offCentre)*condenseRate.oldTime()*rho.oldTime()
      - offCentre*condenseRate*rho
    )
)/rho;

condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
);

