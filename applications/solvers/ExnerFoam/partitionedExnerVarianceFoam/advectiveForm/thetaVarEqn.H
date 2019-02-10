// advective form theta equation
for(label ip = 0; ip < nParts; ip++)
{
    thetaVar.ddt()[ip] =-fvc::div(volFlux[ip], thetaVar[ip], "thetaVar")
                    + thetaVar[ip]*fvc::div(volFlux[ip]);

    // Update theta according to the rate of change
    thetaVar[ip] = thetaVar[ip].oldTime() + dt*
    (
        (1-offCentre)*thetaVar.ddt()[ip].oldTime()
      + offCentre*thetaVar.ddt()[ip]
    );
}

// Mass transfer terms (for nParts = 2 only) block implicit solution
if (nParts != 2)
{
    FatalErrorIn("thetaVarEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

// Update the sum
thetaVar.updateSum();

