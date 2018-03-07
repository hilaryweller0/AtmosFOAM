Info << "Initialising diagnostics file diags.dat\n" << endl;

OFstream diags(args.rootPath() / args.caseName()/"diags.dat");
diags << "#time   maxCo0   maxCo1   minSigma   maxSigma   KE0  KE1 PE0  PE1   IE0  IE1 totalEnergy KEsum PEsum IEsum" << endl;

scalarList maxCoNums(nParts, scalar(0));

const dimensionedScalar Vtot = sum(mesh.V());

scalarList energyInit(nParts, scalar(0));
scalarList KEInit(nParts, scalar(0));
scalarList PEInit(nParts, scalar(0));
scalarList IEInit(nParts, scalar(0));

// Caclulate the initial energy in each partition
for(label ip = 0; ip < partNames.size(); ip++)
{
    volScalarField energyTmp = 0.5*sigmaRho[ip]*magSqr(u[ip]);
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

