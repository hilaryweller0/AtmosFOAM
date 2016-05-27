// Calculate the condensation rate and advect the water vapour and liquid

// Calculate the saturation vapour pressure, es
volScalarField es("es", Pcc*pRef*Foam::exp(-Lv0/Rv*(1/T - 1/T0)));

// Calculation of condensation rate based on Jahn, Knoth, Konig and Vogelsberg
volScalarField condForce = rho*qv - es/(T*Rv);
volScalarField condLim = rho*ql;

condenseRate = 
(
    (condForce - condLim + Foam::sqrt(sqr(condForce) + sqr(condLim)))/dt
);

drhovdt = -fvc::div(U, qv) - condenseRate;
drholdt = -fvc::div(U, ql) + condenseRate;

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

