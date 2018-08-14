fileName invs = args.rootPath() / args.caseName() / "energy.dat";
OFstream os(invs);
os << "#time   normalisedEnergyChange" << endl << "#";

fileName invsPos = args.rootPath() / args.caseName() / "energyPositive.dat";
OFstream osPos(invsPos);

dimensionedScalar energy("energy", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energy0("energy0", dimensionSet(0,3,-2,0,0), scalar(1));
dimensionedScalar energy_old("energy_old", dimensionSet(0,3,-2,0,0), scalar(1e+120));
scalar normalEnergyChange = 0;
scalar deltaNormalEnergyChange = 0;
const dimensionedScalar Vtot = sum(mesh.V());

#include "energy.H"

energy0 = energy;
#include "energy.H"

