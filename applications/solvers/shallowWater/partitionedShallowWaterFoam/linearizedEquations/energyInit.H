fileName invs = args.rootPath() / args.caseName() / "energy.dat";
OFstream os(invs);
os << "#time   normalisedEnergyChange" << endl << "#";

fileName invsPos = args.rootPath() / args.caseName() / "energyPositive.dat";
OFstream osPos(invsPos);

dimensionedScalar energyTot("energyTot", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energyTot0("energyTot0", dimensionSet(0,3,-2,0,0), scalar(1));
dimensionedScalar energy_old("energy_old", dimensionSet(0,3,-2,0,0), scalar(1e+120));
scalar normalEnergyChange = 0;
scalar deltaNormalEnergyChange = 0;
const dimensionedScalar Vtot = sum(mesh.V());

#include "energy.H"

energyTot0 = energyTot;
#include "energy.H"

