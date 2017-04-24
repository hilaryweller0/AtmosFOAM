// Advect densities of all phases of all partitions

//// Update total density
//atmosParts.dRhodt() = -fvc::div(atmosParts.flux());
//atmosParts.rho() = atmosParts.rho().oldTime() + dt*
//(
//    (1-offCentre)*atmosParts.dRhodt().oldTime()
//  + offCentre*atmosParts.dRhodt()
//);

// Loop over all partitions >0 and update partition density
//for(label ip = 1; ip < atmosParts.size(); ip++)
for(label ip = 0; ip < atmosParts.size(); ip++)
{
    partition& parti = atmosParts[ip];
    volScalarField& sigmaRhoi = parti.sigmaRho();
    const surfaceScalarField& fluxi = parti.flux();
    
    // Advect density for this partition (without partition transfer terms)
    parti.dSigmaRhodt() = -fvc::div(fluxi);

    sigmaRhoi = sigmaRhoi.oldTime() + dt*
    (
        (1-offCentre)*parti.dSigmaRhodt().oldTime()
      + offCentre*parti.dSigmaRhodt()
    );
}
//// Update density for partition zero
//partition& part0 = atmosParts[0];
//part0.sigmaRho() = atmosParts.rho();
//for(label ip = 1; ip < atmosParts.size(); ip++)
//{
//    part0.sigmaRho() -= atmosParts[ip].sigmaRho();
//}
//part0.dSigmaRhodt() = -fvc::div(part0.flux());

atmosParts.sumDensity();
atmosParts.dRhodt() = -fvc::div(atmosParts.flux());

// Update and species densities for all partitions
for(label ip = 0; ip < atmosParts.size(); ip++)
{
    partition& parti = atmosParts[ip];
    const volScalarField& sigmaRhoi = parti.sigmaRho();
    const surfaceScalarField& fluxi = parti.flux();
    const volScalarField& sigmai = parti.sigma();

    // Loop over all atmospheric constituents (species) except air (is=0)
    for(label is = 1; is < parti.size(); is++)
    {
        fluidSpecie& specie = parti[is];
        
        // Update the rate of change and advect the gaseous phase of specie
        specie.gas().dRhodt() = -fvc::div
        (
            fluxi, sigmai*specie.gas().rho()/sigmaRhoi
        );

        specie.gas().rho() = 
        (
            sigmai.oldTime()*specie.gas().rho().oldTime() + dt*
            (
                (1-offCentre)*specie.gas().dRhodt().oldTime()
              + offCentre*specie.gas().dRhodt()
            )
        )/sigmai;

        // And the liquid phase
        specie.liquid().dvdt() = -fvc::div
        (
            fluxi, sigmai*specie.liquid().v()/sigmaRhoi
        );

        specie.liquid().v() = 
        (
            sigmai.oldTime()*specie.liquid().v().oldTime() + dt*
            (
                (1-offCentre)*specie.liquid().dvdt().oldTime()
              + offCentre*specie.liquid().dvdt()
            )
        )/sigmai;
    }
    
    // Update air density for this partition from the total
    parti[0].gas().rho() = sigmaRhoi
                         - parti[0].liquid().rho()*parti[0].liquid().v();
    for(label ip = 1; ip < parti.size(); ip++)
    {
        parti[0].gas().rho() -= sigmai*
        (
            parti[ip].gas().rho()
         + parti[ip].liquid().rho()*parti[ip].liquid().v()
        );
    }
    parti[0].gas().rho() /= sigmai;
    
    // Loop over all atmospheric constituents (species) except air (is=0)
    for(label is = 1; is < parti.size(); is++)
    {
        fluidSpecie& specie = parti[is];
        
        // Apply condensation to volatile species (ie water, not air)
        if (specie.pvs0() < specie.gas().p0())
        {
            specie.updateCondensation(parti.T());

            specie.gas().rho()  -= offCentre*specie.condensation();
            specie.liquid().v() += offCentre*specie.condensation()
                                        /specie.liquid().rho();
            specie.gas().dRhodt() -= sigmai*specie.condensation()/dt;
            specie.liquid().dvdt() += sigmai*specie.condensation()/dt
                                        /specie.liquid().rho();
        }
    }
}

