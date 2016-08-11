// Advect all phases
for(label ip = 0; ip < atmos.size(); ip++)
{
    // The gaseous phase
    atmos[ip].gas().dRhodt() = -fvc::div(un, atmos[ip].gas().rho());

    atmos[ip].gas().rho() = atmos[ip].gas().rho().oldTime() + dt*
    (
        (1-offCentre)*atmos[ip].gas().dRhodt().oldTime()
      + offCentre*atmos[ip].gas().dRhodt()
    );

    // The liquid phase
    atmos[ip].liquid().dvdt() = -fvc::div(un, atmos[ip].liquid().v());

    atmos[ip].liquid().v() = atmos[ip].liquid().v().oldTime() + dt*
    (
        (1-offCentre)*atmos[ip].liquid().dvdt().oldTime()
      + offCentre*atmos[ip].liquid().dvdt()
    );
}

// Calculate the condensation for each species
for(label ip = 0; ip < atmos.size(); ip++)
{
    volScalarField condenseRate(atmos[ip].condensation(T)/dt);
    
    atmos[ip].gas().rho()  -= offCentre*dt*condenseRate;
    atmos[ip].liquid().v() += offCentre*dt*condenseRate
                                /atmos[ip].liquid().rho();
    atmos[ip].gas().dRhodt() -= condenseRate;
    atmos[ip].liquid().dvdt() += condenseRate
                                /atmos[ip].liquid().rho();
}
