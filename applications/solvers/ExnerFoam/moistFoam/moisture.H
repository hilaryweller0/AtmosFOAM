// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));
volScalarField rvs = epsilon*es/(p-es);

//// Condensation rate from KJJV15
//// Calculation of condensation rate based on Jahn, Knoth, Konig and Vogelsberg
//volScalarField condForce("condForce", rho*qv - es/(T*Rv));
//volScalarField condLim("condLim", rho*ql);

//condenseRate = 
//(
//    (condForce - condLim + Foam::sqrt(sqr(condForce) + sqr(condLim)))
//   /tauCondense
//);

// Advect mixing ratios
drhovdt = -fvc::div(U, qv);
drholdt = -fvc::div(U, ql);

qv =
(
    qv.oldTime()*rho.oldTime() + dt*
    (
        offCentre*drhovdt + (1-offCentre)*drhovdt.oldTime()
    )
)/rho;

ql =
(
    ql.oldTime()*rho.oldTime() + dt*
    (
        offCentre*drholdt + (1-offCentre)*drholdt.oldTime()
    )
)/rho;

// Update related mixing ratios
rv = qv/qd;
rl = ql/qd;

// Calculate and apply the condensation
condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
);

// Apply condensation
drhovdt -= rho*condenseRate;
drholdt += rho*condenseRate;

qv -= dt*offCentre*condenseRate;
ql += dt*offCentre*condenseRate;

// Update related mixing ratios
qd = 1 - qv - ql;
rv = qv/qd;
rl = ql/qd;
