// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", pSat(T, pSat0, Lv0, Rv, T0));

//// Condensation rate from KJJV15
//// Calculation of condensation rate based on Jahn, Knoth, Konig and Vogelsberg
//volScalarField condForce("condForce", rho*qv - es/(T*Rv));
//volScalarField condLim("condLim", rho*ql);

//condenseRate = 
//(
//    (condForce - condLim + Foam::sqrt(sqr(condForce) + sqr(condLim)))
//   /tauCondense
//);

// Update mixing ratios without condensation
drhovdt = -fvc::div(U, qv); // - condenseRate;
drholdt = -fvc::div(U, ql); // + condenseRate;

// Update mass mixing ratios without splitting
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
qd = 1 - qv - ql;
rv = qv/qd;
rl = ql/qd;

// Condensation rate from BF02
volScalarField rvs = epsilon*es/(p-es);
condenseRate = max
(
    (rv - rvs)/(dt*(1 + sqr(Lv)*rvs/(Cp*Rv*sqr(T)))),
    -rl/dt
)*rho;

// Apply condensation
rv -= dt*condenseRate/rho;
rl += dt*condenseRate/rho;
qv = rv/(1+rv+rl);
ql = rl/(1+rv+rl);
qd = 1 - qv - ql;

