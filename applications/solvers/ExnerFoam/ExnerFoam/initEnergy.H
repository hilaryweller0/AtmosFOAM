Info << "Initialising total energy file energy.dat\n" << endl;

OFstream es(args.rootPath() / args.caseName() / "energy.dat");
es << "#time   KE       internal potential normalisedEnergy  max w max theta"
   << endl;

dimensionedScalar Vtot = sum(mesh.V());

volScalarField KE
(
    IOobject("KE", runTime.timeName(), mesh),
    0.25*fvc::surfaceSum(sqr(phi)/(mesh.magSf()*mesh.deltaCoeffs()*rhof))/dimensionedScalar("V", dimVol, scalar(1))
);
KE.internalField() /= mesh.V();
volScalarField internalEnergy
(
    IOobject("internalEnergy", runTime.timeName(), mesh),
    rho*theta*Exner*Cv
);
volScalarField potentialEnergy
(
    IOobject("potentialEnergy", runTime.timeName(), mesh),
    -rho*(g & mesh.C())
);

const dimensionedScalar KEInitTot = fvc::domainIntegrate(KE)/Vtot;
const dimensionedScalar internalInitTot = fvc::domainIntegrate(internalEnergy)/Vtot;
const dimensionedScalar potentialInitTot = fvc::domainIntegrate(potentialEnergy)/Vtot;

const dimensionedScalar energyInitTot = KEInitTot + internalInitTot + potentialInitTot;

