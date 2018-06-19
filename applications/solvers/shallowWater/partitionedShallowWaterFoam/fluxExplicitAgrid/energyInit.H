fileName invs = args.rootPath() / args.caseName() / "energy.dat";
OFstream os(invs);
os << "#time   KE    PE      normalisedEnergyChange" << endl << "#";

dimensionedScalar KE("KE", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar PE("PE", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energy("energy", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar energy0("energy0", dimensionSet(0,3,-2,0,0), scalar(1));
scalar normalEnergyChange = 0;
const dimensionedScalar Vtot = sum(mesh.V());

#include "energy.H"

energy0 = energy;
#include "energy.H"

