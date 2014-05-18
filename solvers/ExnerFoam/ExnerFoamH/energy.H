{
    KE = 0.25*fvc::surfaceSum(U*V/rhof)
        /dimensionedScalar("V", dimVol, scalar(1));
    KE.internalField() /= mesh.V();
    internalEnergy = rho*theta*Exner*Cv;
    potentialEnergy = -rho*(g & mesh.C());

    dimensionedScalar KETot = fvc::domainIntegrate(KE)/Vtot;
    dimensionedScalar internalTot = fvc::domainIntegrate(internalEnergy)/Vtot;
    dimensionedScalar potentialTot = fvc::domainIntegrate(potentialEnergy)/Vtot;
    dimensionedScalar energyDiffNorm
         = (KETot + internalTot + potentialTot - energyInitTot)/energyInitTot;
    
    es << runTime.timeName() << "  "
       << ((KETot-KEInitTot)/energyInitTot).value() << "  "
       << ((internalTot-internalInitTot)/energyInitTot).value() << "  "
       << ((potentialTot-potentialInitTot)/energyInitTot).value() << "  "
       << energyDiffNorm.value() << "  "
       << max(Uf.component(vector::Z)).value() << "  "
       << max(theta).value() << endl;
}

