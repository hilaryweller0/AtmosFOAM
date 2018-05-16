// Potential energy
energyField =  0.5*g*sqr(h);

PE1 = fvc::domainIntegrate(energyField)/Vtot;

// Kinetic energy
energyField *= 0;
for(label ip = 0; ip < nParts; ip++)
{
    energyField += 0.5*sigmah[ip]*magSqr(u[ip]);
}

KE1 = fvc::domainIntegrate(energyField)/Vtot;

E1 = KE1 + PE1;

F1 *= 0;
for(label ip = 0; ip < nParts; ip++)
{
    F1 += fvc::interpolate(sigmah[ip])*volFlux[ip];
    //F1 += flux[ip];
}
