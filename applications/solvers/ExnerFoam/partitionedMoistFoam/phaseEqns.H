// Advect densities of all phases of all partitions

// Loop over all partitions and update partition density and species densities
for(label ip = 0; ip < atmosParts.size(); ip++)
{
    partition& parti = atmosParts[ip];
    volScalarField& rhoi = parti.rho();
    surfaceScalarField& fluxi = parti.flux();
    
    // Advect density for this partition
    // (without partition transfer terms)
    parti.dRhodt() = -fvc::div(fluxi);

    rhoi = rhoi.oldTime() + dt*
    (
        (1-offCentre)*parti.dRhodt().oldTime()
      + offCentre*parti.dRhodt()
    );

    // Loop over all atmospheric constituents (species) except air (is=0)
    for(label is = 1; is < parti.size(); is++)
    {
        fluidSpecie& specie = parti[is];
        
        // Update the rate of change and advect the gaseous phase of specie
        specie.gas().dRhodt() = -fvc::div(fluxi, specie.gas().rho()/rhoi);

        specie.gas().rho() = specie.gas().rho().oldTime() + dt*
        (
            (1-offCentre)*specie.gas().dRhodt().oldTime()
          + offCentre*specie.gas().dRhodt()
        );

        // And the liquid phase
        specie.liquid().dvdt() = -fvc::div(fluxi, specie.liquid().v()/rhoi);

        specie.liquid().v() = specie.liquid().v().oldTime() + dt*
        (
            (1-offCentre)*specie.liquid().dvdt().oldTime()
          + offCentre*specie.liquid().dvdt()
        );
        
        // Apply condensation to volatile species (ie water, not air)
        if (specie.pvs0() < specie.gas().p0())
        {
            //specie.updateCondensation(T);
        
            specie.gas().rho()  -= offCentre*specie.condensation();
            specie.liquid().v() += offCentre*specie.condensation()
                                        /specie.liquid().rho();
            specie.gas().dRhodt() -= specie.condensation()/dt;
            specie.liquid().dvdt() += specie.condensation()/dt
                                        /specie.liquid().rho();
        }
    }

    // Update air density for this partition
    parti[0].gas().rho() += rhoi - parti.sumDensity();
}

// Update total density
atmosParts.sumDensity();

