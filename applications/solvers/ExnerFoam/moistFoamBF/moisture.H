// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));

// Condensation rate from BF02
volScalarField rvs = epsilon*es/(p-es);
condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
)*rho;

// Update rv and rl
drhovdt = -fvc::div(U, rv) - condenseRate;
drholdt = -fvc::div(U, rl) + condenseRate;

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

