{
    // Potential energy
    energyField =  0.5*g*sqr(h);
    
    PE2 = fvc::domainIntegrate(energyField)/Vtot;
    
    // Kinetic energy
    energyField *= 0;
    for(label ip = 0; ip < nParts; ip++)
    {
        energyField += 0.5*sigmah[ip]*magSqr(u[ip]);
        delta_ke += 0.5*fvc::interpolate(sigmah[ip])*volFlux[ip]*volFlux[ip];
    }

    KE2 = fvc::domainIntegrate(energyField)/Vtot;

    E2 = KE2 + PE2;
    
    deltaKE = KE2 - KE1;
    
    F2 *= 0;
    for(label ip = 0; ip < nParts; ip++)
    {
        F2 += fvc::interpolate(sigmah[ip])*volFlux[ip];
        //F2 += flux[ip];
    }
    
    deltaF = mag(fvc::domainIntegrate(fvc::reconstruct(F2-F1))) / Vtot;
    
    Info << "Change in momentum min/max: " << min(F2-F1).value() << " " << max(F2-F1).value() << endl;
    Info << "Change in kinetic energy: " << deltaKE.value() << endl;
    

    energyChangeStream << runTime.timeName() << " " 
                       << KE2.value() - KE1.value() << " " 
                       << PE2.value() - PE1.value() << " "
                       << deltaF.value() << endl;
}
