Info << "Initialising diagnostics file diags.dat\n" << endl;

OFstream diags(args.rootPath() / args.caseName()/"diags.dat");
diags << "#time   maxCo0   maxCo1   minSigma   maxSigma   energy0   energy1   energySum  energyMean" << endl;

scalarList maxCoNums(nParts, scalar(0));

partitionedVolScalarField energy
(
    "energy",
    partNames,
    0.5*sigmaRho[0]*magSqr(u[0])
  + sigmaRho[0]*theta[0]*Exner*Cv
  - sigmaRho[0]*(g & mesh.C())
);

for(label ip = 1; ip < partNames.size(); ip++)
{
    energy[ip] = 0.5*sigmaRho[0]*magSqr(u[ip])
               + sigmaRho[ip]*theta[ip]*Exner*Cv
               - sigmaRho[ip]*(g & mesh.C());
}

rho.updateSum();
volScalarField energyMean = 0.5*rho.sum()*magSqr(u.updateSum()/rho.sum())
                          + theta.updateSum()*Exner*Cv
                          - rho.sum()*(g & mesh.C());

const dimensionedScalar Vtot = sum(mesh.V());

scalarList energyInit(nParts, scalar(0));

for(label ip = 0; ip < partNames.size(); ip++)
{
    energyInit[ip] = (fvc::domainIntegrate(energy[ip])/Vtot).value();
}
const scalar energyInitTot = sum(energyInit);
const scalar energyMeanInit = (fvc::domainIntegrate(energyMean)/Vtot).value();

