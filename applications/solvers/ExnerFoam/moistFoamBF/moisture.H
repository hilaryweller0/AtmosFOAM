{
// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));
volScalarField rvs = epsilon*es/(p-es);

// Advect mixing ratios
drhovdt = -fvc::div(U, rv);
drholdt = -fvc::div(U, rl);

rv =
(
    rv.oldTime()*rho.oldTime() + dt*
    (
        offCentre*drhovdt + (1-offCentre)*drhovdt.oldTime()
    )
)/rho;

rl =
(
    rl.oldTime()*rho.oldTime() + dt*
    (
        offCentre*drholdt + (1-offCentre)*drholdt.oldTime()
    )
)/rho;

// Calculate and apply the condensation
condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
);

// Apply condensation
drhovdt -= rho*condenseRate;
drholdt += rho*condenseRate;

rv -= dt*offCentre*condenseRate;
rl += dt*offCentre*condenseRate;
}
