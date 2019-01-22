{
// Energy calculations

scalarList energy(nParts, scalar(0));
scalarList KE(nParts, scalar(0));
scalarList PE(nParts, scalar(0));
scalarList IE(nParts, scalar(0));
totalFlux *= 0;
totalMomentum *= 0;
energyField *= 0;

// Caclulate the initial energy in each partition
for(label ip = 0; ip < partNames.size(); ip++)
{
    volScalarField energyTmp = 0.5*sigmaRho[ip]*magSqr(u[ip]);
    energyField += energyTmp;
    KE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyTmp = - sigmaRho[ip]*(g & mesh.C());
    energyField += energyTmp;
    PE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyTmp = sigmaRho[ip]*theta[ip]*sigmaRho.sum()/Psi*Cv;
    energyField += energyTmp;
    IE[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energy[ip] = KE[ip] + PE[ip] + IE[ip];
    
    totalFlux += linearInterpolate(sigmaRho[ip])*volFlux[ip];
    totalMomentum += sigmaRho[ip]*u[ip];
}

scalar KEsum = sum(KE);
scalar PEsum = sum(PE);
scalar IEsum = sum(IE);
scalar energyTot = sum(energy);

diags << runTime.timeName() << "  "
      << maxCoNums[0] << "  " << maxCoNums[1] << "  "
      << min(sigma[1]).value() << "  "
      << max(sigma[1]).value() << "  "
      << (KE[0] - KEInit[0])/energyInitTot << "  "
      << (KE[1] - KEInit[1])/energyInitTot << "  "
      << (PE[0] - PEInit[0])/energyInitTot << "  "
      << (PE[1] - PEInit[1])/energyInitTot << "  "
      << (IE[0] - IEInit[0])/energyInitTot << "  "
      << (IE[1] - IEInit[1])/energyInitTot << "  "
      << (energyTot - energyInitTot)/energyInitTot << "  "
      << (KEsum - KEsumInit)/energyInitTot << "  "
      << (PEsum - PEsumInit)/energyInitTot << "  "
      << (IEsum - IEsumInit)/energyInitTot << endl;
      
scalar fluxChange = sum
(
    (totalFlux - totalFluxPreTransfer) /
    (totalFluxPreTransfer + dimensionedScalar("",totalFlux.dimensions(),scalar(1e-16)))
).value();
scalar momentumChange = 
(
    mag
    (
        fvc::domainIntegrate
        (
            totalMomentum - totalMomentumPreTransfer
        )
    ).value() / 
    (mag
    (
        fvc::domainIntegrate(totalMomentum)
    ).value() + 1e-16)
);

Info    << "Transfer energy change: " 
        << (energyTot - energyTotPreTransfer) << endl;
Info    << "Momentum change: " 
        << fluxChange << endl;

transferDiags 
      << runTime.timeName() << "  "
      << (energyTot - energyTotPreTransfer) << "  "
      << (KEsum - KEsumPreTransfer) << "  "
      << (PEsum - PEsumPreTransfer) << "  "
      << (IEsum - IEsumPreTransfer) << " "
      << momentumChange << " "
      << fluxChange << " "
      << min(energyField - energyFieldPreTransfer).value() << " "
      << max(energyField - energyFieldPreTransfer).value()
      << endl;
}
