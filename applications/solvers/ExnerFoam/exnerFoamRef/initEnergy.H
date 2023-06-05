Info << "Initialising total energy file energy.dat\n" << endl;

fileName OFname = Pstream::parRun() ?
                args.rootPath() / args.caseName() / "../energy.dat" :
                args.rootPath() / args.caseName() / "energy.dat";
OFstream es(OFname);

dimensionedScalar Vtot = sum(mesh.V());

volScalarField KE
(
    IOobject("KE", runTime.timeName(), mesh),
    0.5*rho*magSqr(U)
);
volScalarField internalEnergy
(
    IOobject("internalEnergy", runTime.timeName(), mesh),
    rho*thermo.T()*Cv
);
volScalarField potentialEnergy
(
    IOobject("potentialEnergy", runTime.timeName(), mesh),
    -rho*(g & mesh.C())
);

Info << "KE dimensions = " << KE.dimensions() << nl
     << "internalEnergy dimensions = " << internalEnergy.dimensions() << nl
     << "potentialEnergy dimensions = " << potentialEnergy.dimensions() << endl;

const dimensionedScalar KEInitTot = fvc::domainIntegrate(KE)/Vtot;
const dimensionedScalar internalInitTot 
      = fvc::domainIntegrate(internalEnergy)/Vtot;
const dimensionedScalar potentialInitTot 
      = fvc::domainIntegrate(potentialEnergy)/Vtot;
      
Info << "Vtot = " << Vtot << "\nKEInitTot = " << KEInitTot
     << "\ninternalInitTot = " << internalInitTot << "\npotentialInitTot = "
     << potentialInitTot << endl;
     
const dimensionedScalar energyInitTot 
      = KEInitTot + internalInitTot + potentialInitTot;

if (Pstream::master())
{
    es << "#time   KE       internal potential normalisedEnergy  max w"
       << endl;
}

