// Advect multi-phase phases
for(label ip = 1; ip < atmos.size(); ip++)
{
    // The gaseous phase
    atmos[ip].gas().dRhodt() = -fvc::div(U, atmos[ip].gas().rho()/rho);

    atmos[ip].gas().rho() = atmos[ip].gas().rho().oldTime() + dt*
    (
        (1-offCentre)*atmos[ip].gas().dRhodt().oldTime()
      + offCentre*atmos[ip].gas().dRhodt()
    );

    // The liquid phase
    atmos[ip].liquid().dvdt() = -fvc::div(U, atmos[ip].liquid().v()/rho);

    atmos[ip].liquid().v() = atmos[ip].liquid().v().oldTime() + dt*
    (
        (1-offCentre)*atmos[ip].liquid().dvdt().oldTime()
      + offCentre*atmos[ip].liquid().dvdt()
    );
}

// Update the air density
air.rho() += rho - atmos.sumDensity();

// Calculate and apply the condensation for each species
for(label ip = 0; ip < atmos.size(); ip++)
{
    // Only apply condensation if the vapour pressure is close to ambient
    if (atmos[ip].pvs0() < atmos[ip].gas().p0())
    {
        atmos[ip].updateCondensation(T);
    
        atmos[ip].gas().rho()  -= offCentre*atmos[ip].condensation();
        atmos[ip].liquid().v() += offCentre*atmos[ip].condensation()
                                    /atmos[ip].liquid().rho();
        atmos[ip].gas().dRhodt() -= atmos[ip].condensation()/dt;
        atmos[ip].liquid().dvdt() += atmos[ip].condensation()/dt
                                    /atmos[ip].liquid().rho();
    }
}

//Info << "After advection air.rho goes from "
//     << min(air.rho().internalField()).value() << " to "
//     << max(air.rho().internalField()).value() << nl
//     << "Water vapour density goes from "
//     << min(water.gas().rho().internalField()).value() << " to "
//     << max(water.gas().rho().internalField()).value() << nl
//     << "Water liquid fraction goes from "
//     << min(water.liquid().v().internalField()).value() << " to "
//     << max(water.liquid().v().internalField()).value() << endl;

