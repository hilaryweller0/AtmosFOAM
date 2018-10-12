{
// Energy calculations

scalarList energy(nParts, scalar(0));
scalarList KE(nParts, scalar(0));
scalarList PE(nParts, scalar(0));
scalarList IE(nParts, scalar(0));

// Caclulate the initial energy in each partition
for(label ip = 0; ip < partNames.size(); ip++)
{
    volScalarField energyTmp = 0.5*sigmaRho[ip]*magSqr(u[ip]);
    KE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyTmp = - sigmaRho[ip]*(g & mesh.C());
    PE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyTmp = sigmaRho[ip]*theta[ip]*sigmaRho.sum()/Psi*Cv;
    IE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energy[ip] = KE[ip] + PE[ip] + IE[ip];
}

KEsumPreTransfer = sum(KE);
PEsumPreTransfer = sum(PE);
IEsumPreTransfer = sum(IE);
energyTotPreTransfer = sum(energy);
}
