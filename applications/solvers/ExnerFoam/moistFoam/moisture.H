// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", Pcc*pRef*Foam::exp(TccScale*(T - T0)/(T - Tcc)));

// Calculate the saturation mixing ratio
volScalarField rvs = epsilon*(es/(p-es));
volScalarField Lv("Lv", Lv0 - (Cpl - Cpv)*(T - T0));


rv = rv.oldTime() - dt*(1-offCentre)*condenseRate.oldTime();
rl = rl.oldTime() + dt*(1-offCentre)*condenseRate.oldTime();

rv =
(
    rv.oldTime()*rho.oldTime() - dt*
    (
        (1-offCentre)*divUrv.oldTime()
      + offCentre*fvc::div(U,rv)
    )
)/rho;

rl =
(
    rl.oldTime()*rho.oldTime() - dt*
    (
        (1-offCentre)*divUrl.oldTime()
      + offCentre*fvc::div(U,rl)
    )
)/rho;

condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
);

rv -= dt*offCentre*condenseRate;
rl += dt*offCentre*condenseRate;

