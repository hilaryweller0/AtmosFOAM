fileName energyChangeDat = args.rootPath() / args.caseName() / "energyChangeTransfers.dat";
OFstream energyChangeStream(energyChangeDat);
energyChangeStream << "#time   delta_KE    delta_PE    delta_F" << endl;

volScalarField energyField("energyField", 0.5*g*sqr(h));
dimensionedScalar E1("E1", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar KE1("KE1", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar PE1("PE1", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar E2("E2", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar KE2("KE2", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar PE2("PE2", dimensionSet(0,3,-2,0,0), scalar(0));
dimensionedScalar deltaF("deltaF", dimensionSet(0,2,-1,0,0), scalar(0));
surfaceScalarField F1("F1", fvc::interpolate(sigmah[0])*volFlux[0]);
surfaceScalarField F2("F2", fvc::interpolate(sigmah[0])*volFlux[0]);
//const dimensionedScalar Vtot = sum(mesh.V());

