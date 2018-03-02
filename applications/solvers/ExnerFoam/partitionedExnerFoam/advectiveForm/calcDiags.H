// Energy calculations

for(label ip = 0; ip < partNames.size(); ip++)
{
    energy[ip] = 0.5*sigmaRho[0]*magSqr(u[ip])
               + sigmaRho[ip]*theta[ip]*Exner*Cv
               - sigmaRho[ip]*(g & mesh.C());
}

energyMean = 0.5*rho.sum()*magSqr(u.updateSum()/rho.sum())
           + theta.updateSum()*Exner*Cv
           - rho.sum()*(g & mesh.C());

scalar energyMeanSum = (fvc::domainIntegrate(energyMean)/Vtot).value();


scalarList energySum(nParts, scalar(0));

for(label ip = 0; ip < partNames.size(); ip++)
{
    energySum[ip] = (fvc::domainIntegrate(energy[ip])/Vtot).value();
}

const scalar energyTot = sum(energySum);

diags << runTime.timeName() << "  "
      << maxCoNums[0] << "  " << maxCoNums[1] << "  "
      << min(sigma[1]).value() << "  "
      << max(sigma[1]).value() << "  "
      << (energySum[0]-energyInit[0])/energyInitTot << "  "
      << (energySum[1]-energyInit[1])/energyInitTot << "  "
      << (energyTot-energyInitTot)/energyInitTot << "  "
      << (energyMeanSum - energyMeanInit)/energyMeanInit << endl;

