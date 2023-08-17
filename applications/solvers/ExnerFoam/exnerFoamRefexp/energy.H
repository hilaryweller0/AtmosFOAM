{
    KE = 0.5*rho*magSqr(U);
    internalEnergy = rho*thermo.T()*Cv;
    potentialEnergy = -rho*(g & mesh.C());

    dimensionedScalar KETot = fvc::domainIntegrate(KE)/Vtot;
    dimensionedScalar internalTot = fvc::domainIntegrate(internalEnergy)/Vtot;
    dimensionedScalar potentialTot = fvc::domainIntegrate(potentialEnergy)/Vtot;
    dimensionedScalar energyDiffNorm
         = (KETot + internalTot + potentialTot - energyInitTot)/energyInitTot;
    dimensionedScalar massTot = fvc::domainIntegrate(thermo.p()/(R*thermo.T()));
    
    scalar maxw = max(U.component(vector::Z)).value();
    
    if (Pstream::master())
    {
        es << runTime.timeName() << "  "
           << ((KETot-KEInitTot)/energyInitTot).value() << "  "
           << ((internalTot-internalInitTot)/energyInitTot).value() << "  "
           << ((potentialTot-potentialInitTot)/energyInitTot).value() << "  "
           << energyDiffNorm.value() << " "
           << maxw << " " << ((massTot - massTot0)/massTot0).value() << endl;
    }
}

