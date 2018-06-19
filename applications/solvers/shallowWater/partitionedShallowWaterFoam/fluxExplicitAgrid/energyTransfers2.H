{
    // Potential energy
    energyField =  0.5*g*sqr(h);
    
    PE2 = fvc::domainIntegrate(energyField)/Vtot;
    
    // Kinetic energy
    energyField *= 0;
    for(label ip = 0; ip < nParts; ip++)
    {
        energyField += 0.5*sigmah[ip]*magSqr(u[ip]);
        delta_ke_centre += 0.5*sigmah[ip]*magSqr(u[ip]);
    }

    KE2 = fvc::domainIntegrate(energyField)/Vtot;

    E2 = KE2 + PE2;
    
    deltaKE = KE2 - KE1;
    
    hu2 *= 0;
    for(label ip = 0; ip < nParts; ip++)
    {
        hu2 += hu[ip];
    }
    
    deltaF = mag(fvc::domainIntegrate(mag(hu2-hu1))) / Vtot;
    
    Info << "Change in momentum min/max: " << min(mag(hu2-hu1)).value() << " " << max(mag(hu2-hu1)).value() << endl;
    Info << "Change in kinetic energy: " << deltaKE.value() << endl;
    Info << "Min/Max KE change Cell Centre: " << min(delta_ke_centre).value() << " " << max(delta_ke_centre).value() << endl;
    

    energyChangeStream << runTime.timeName() << " " 
                       << KE2.value() - KE1.value() << " " 
                       << PE2.value() - PE1.value() << " "
                       << deltaF.value() << endl;
}
