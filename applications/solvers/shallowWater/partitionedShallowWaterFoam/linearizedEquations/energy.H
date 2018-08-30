{
    // Potential energy
    energy *= 0;
    energy += 0.5*gravity*sqr(sigmah[0] + sigmah[1]);
    // Kinetic energy
    for(label ip = 0; ip < nParts; ip++)
    {
        energy += 0.5*SIGMAH[ip]*magSqr(u[ip]);
        energy +=     sigmah[ip]*(U[ip] & u[ip]);
    }

    energyTot = fvc::domainIntegrate(energy)/Vtot;

    normalEnergyChange = ((energyTot - energyTot0)/energyTot0).value();
    deltaNormalEnergyChange = ((energyTot - energy_old)/energyTot0).value();
    
    os << runTime.timeName() << " " << normalEnergyChange << endl;
       
    if ( deltaNormalEnergyChange > 0. )
    {
        osPos << runTime.timeName() << " " << deltaNormalEnergyChange << endl;
    }
    energy_old = energyTot;
}
