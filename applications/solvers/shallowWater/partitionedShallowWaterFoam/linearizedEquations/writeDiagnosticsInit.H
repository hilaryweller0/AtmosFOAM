fileName file1 = args.rootPath() / args.caseName() / "minmax.dat";
OFstream minmax(file1);
minmax << "#time   min h    max h   min sigma   max sigma " << endl;

fileName file2 = args.rootPath() / args.caseName() / "total.dat";
OFstream total(file2);
total << "#time   h    sigma1 " << endl;

#include "writeDiagnostics.H"
