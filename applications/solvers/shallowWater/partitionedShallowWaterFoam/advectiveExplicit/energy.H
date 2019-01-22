{
    // Potential energy
    volScalarField energyField("energyField", 0.5*gravity*sqr(h));
    PE = fvc::domainIntegrate(energyField)/Vtot;
    
    energy = energyField;
    
    // Kinetic energy
    energyField *= 0;
    for(label ip = 0; ip < nParts; ip++)
    {
        energyField += 0.5*sigmah[ip]*magSqr(u[ip]);
    }
    KE = fvc::domainIntegrate(energyField)/Vtot;
    energyTot = KE + PE;

    energy += energyField;

    normalEnergyChange = ((energyTot - energyTot0)/energyTot0).value();
    deltaNormalEnergyChange = ((energyTot - energy_old)/energyTot0).value();
    
    os << runTime.timeName() << " " << KE.value() << " " << PE.value() << " "
       << normalEnergyChange << endl;
       
    if ( deltaNormalEnergyChange > 0. )
    {
        osPos << runTime.timeName() << " " << deltaNormalEnergyChange << endl;
    }
    energy_old = energyTot;
}
