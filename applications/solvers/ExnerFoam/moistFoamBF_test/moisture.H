// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", Pcc*pRef*Foam::exp(-Lv0/Rv*(1/T - 1/T0)));

// Calculate the saturation mixing ratio
volScalarField rvs = epsilon*es/(p-es);
Lv = Lv0 - (Cpl - Cpv)*(T - T0);

rvSource = -fvc::div(U, rv) + rv*divu - condenseRate;
rlSource = -fvc::div(U, rl) + rl*divu + condenseRate;

// Update rv and rl without splitting
rv = rv.oldTime()
   + dt*(1-offCentre)*rvSource.oldTime()
   + dt*offCentre*rvSource;

rl = rl.oldTime()
   + dt*(1-offCentre)*rlSource.oldTime()
   + dt*offCentre*rlSource;

condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
);

