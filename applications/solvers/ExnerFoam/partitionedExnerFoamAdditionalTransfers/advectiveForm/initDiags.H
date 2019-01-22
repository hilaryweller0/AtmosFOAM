Info << "Initialising diagnostics file diags.dat\n" << endl;

OFstream diags(args.rootPath() / args.caseName()/"diags.dat");
diags << "#time   maxCo0   maxCo1   minSigma   maxSigma   KE0  KE1 PE0  PE1   IE0  IE1 totalEnergy KEsum PEsum IEsum" << endl;

scalarList maxCoNums(nParts, scalar(0));

const dimensionedScalar Vtot = sum(mesh.V());

scalarList energyInit(nParts, scalar(0));
scalarList KEInit(nParts, scalar(0));
scalarList PEInit(nParts, scalar(0));
scalarList IEInit(nParts, scalar(0));

// Read in the backgroud potential temperature and calculate its
// internal energy
volScalarField thetaBack
(
    IOobject
    (
        "theta_init",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);
volScalarField energyTmp = sigmaRho.sum()*thetaBack*Exner*Cv;
const scalar IEBack = (fvc::domainIntegrate(energyTmp)/Vtot).value();

// Caclulate the initial energy in each partition
for(label ip = 0; ip < partNames.size(); ip++)
{
    energyTmp = 0.5*sigmaRho[ip]*magSqr(u[ip]);
    KEInit[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    energyTmp = - sigmaRho[ip]*(g & mesh.C());
    PEInit[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    energyTmp = sigmaRho[ip]*theta[ip]*Exner*Cv;
    IEInit[ip] = (fvc::domainIntegrate(energyTmp)/Vtot).value();
    
    energyInit[ip] = KEInit[ip] + PEInit[ip] + IEInit[ip];
}

const scalar KEsumInit = sum(KEInit);
const scalar PEsumInit = sum(PEInit);
const scalar IEsumInit = sum(IEInit);
const scalar energyInitTot = sum(energyInit);
const scalar availableEnergy = IEsumInit - IEBack;
Info << "Available energy = " << availableEnergy << endl;
Info << "energyInitTot = " << energyInitTot << endl;
Info << "Initial KE = " << KEsumInit << endl;
Info << "Initial PE = " << PEsumInit << endl;
Info << "Initial IE = " << IEsumInit << endl;

Info << "Initialising diagnostics file transferDiags.dat\n" << endl;

OFstream transferDiags(args.rootPath() / args.caseName()/"transferDiags.dat");
transferDiags << "#time   deltaTotalEnergy deltaKEsum deltaPEsum deltaIEsum"
              << " deltaMomentum deltaFlux minEnergyChange maxEnergyChange" 
              << endl;
              
scalar KEsumPreTransfer = 0;
scalar PEsumPreTransfer = 0;
scalar IEsumPreTransfer = 0;
scalar energyTotPreTransfer = 0;

surfaceScalarField totalFluxPreTransfer = 0*flux[0];
surfaceScalarField totalFlux = 0*flux[0];
volVectorField totalMomentumPreTransfer = 0*sigmaRho[0]*u[0];
volVectorField totalMomentum = 0*sigmaRho[0]*u[0];
volScalarField energyFieldPreTransfer = 0*sigmaRho[0]*magSqr(u[0]);
volScalarField energyField = 0*sigmaRho[0]*magSqr(u[0]);

