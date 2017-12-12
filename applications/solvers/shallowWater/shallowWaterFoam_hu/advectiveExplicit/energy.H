{
    // Potential energy
    volScalarField energyField("energyField", 0.5*g*sqr(h));
    
    PE = fvc::domainIntegrate(energyField)/Vtot;
    
    // Kinetic energy
    energyField = 0.5*h*magSqr(U);
    KE = fvc::domainIntegrate(energyField)/Vtot;

    energy = KE + PE;

    normalEnergyChange = ((energy - energy0)/energy0).value();
    
    os << runTime.timeName() << " " << KE.value() << " " << PE.value() << " "
       << normalEnergyChange << endl;
}
