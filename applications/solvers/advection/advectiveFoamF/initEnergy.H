Info << "Initialising total energy file energy.dat\n" << endl;

OFstream es(args.rootPath() / args.caseName() / "energy.dat");
es << "#time normalisedSum" << endl;

scalar totalA = 0;
scalar totalInitialSum = 0;

forAll(Tf, facei)
{
    const scalar value = Tf[facei];
    const scalar area = mesh.delta()()[facei] & mesh.Sf()[facei];
    totalA += area;
    totalInitialSum += value*area;
}

totalInitialSum /= totalA;
