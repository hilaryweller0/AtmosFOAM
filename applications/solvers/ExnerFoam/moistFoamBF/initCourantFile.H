OFstream courantFile(args.rootPath() / args.caseName() / "courant.dat");
courantFile << "#time meanCourantNo maxCourantNo Ndt" << endl;
