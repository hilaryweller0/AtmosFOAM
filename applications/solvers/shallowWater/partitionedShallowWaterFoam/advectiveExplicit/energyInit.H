fileName invs = args.rootPath() / args.caseName() / "energy.dat";
OFstream os(invs);
os << "#time   KE    PE      normalisedEnergyChange" << endl << "#";

fileName invsPos = args.rootPath() / args.caseName() / "energyPositive.dat";
OFstream osPos(invsPos);

dimensionedScalar KE("KE", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar PE("PE", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energyTot("energy", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energyTot0("energy0", dimensionSet(0,3,-2,0,0), scalar(1));
dimensionedScalar energy_old("energy_old", dimensionSet(0,3,-2,0,0), scalar(1e+120));
scalar normalEnergyChange = 0;
scalar deltaNormalEnergyChange = 0;
const dimensionedScalar Vtot = sum(mesh.V());

#include "energy.H"

energyTot0 = energyTot;
#include "energy.H"

