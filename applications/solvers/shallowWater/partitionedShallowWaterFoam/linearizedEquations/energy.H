{
    // Potential energy
    volScalarField energyField("energyField", 0.5*gravity*sqr(h));
    
    energyField *= 0;
    energyField += 0.5*gravity*sqr(sigmah[0] + sigmah[1]);
    for(label ip = 0; ip < nParts; ip++)
    {
        energyField += 0.5*SIGMAH[ip]*magSqr(u[ip]);
        energyField +=     sigmah[ip]*(U[ip] & u[ip]);
    }

    energy = fvc::domainIntegrate(energyField)/Vtot;

    normalEnergyChange = ((energy - energy0)/energy0).value();
    deltaNormalEnergyChange = ((energy - energy_old)/energy0).value();
    
    os << runTime.timeName() << " " << normalEnergyChange << endl;
       
    if ( deltaNormalEnergyChange > 0. )
    {
        osPos << runTime.timeName() << " " << deltaNormalEnergyChange << endl;
    }
    energy_old = energy;
}
