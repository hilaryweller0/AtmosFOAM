Info << "Initialising total energy file energy.dat\n" << endl;

OFstream es(args.rootPath() / args.caseName() / "energy.dat");
es << "#time normalisedSum" << endl;

const dimensionedScalar totalV = sum(mesh.V());
const dimensionedScalar totalInitialSum = fvc::domainIntegrate(T)/totalV;
