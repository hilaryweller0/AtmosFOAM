//Minima and maxima of fields.
fileName file1 = args.rootPath() / args.caseName() / "minmax_prognostic.dat";
OFstream minmax_prognostic(file1);
minmax_prognostic << "#time   min rl    max rl   min rv   max rv   min rt   max rt " << endl;

fileName file2 = args.rootPath() / args.caseName() / "minmax_diagnostic.dat";
OFstream minmax_diagnostic(file2);
minmax_diagnostic << "#time   min rl    max rl   min rv   max rv   min rt   max rt " << endl;

fileName file3 = args.rootPath() / args.caseName() / "minmax_analytic.dat";
OFstream minmax_analytic(file3);
minmax_analytic << "#time   min rl    max rl   min rv   max rv   min rt   max rt " << endl;


//Conservation properties of fields.
fileName file4 = args.rootPath() / args.caseName() / "total_prognostic.dat";
OFstream total_prognostic(file4);
total_prognostic << "#time   rl    rv    rt" << endl;

fileName file5 = args.rootPath() / args.caseName() / "total_diagnostic.dat";
OFstream total_diagnostic(file5);
total_diagnostic << "#time   rl    rv    rt" << endl;

fileName file6 = args.rootPath() / args.caseName() / "total_analytic.dat";
OFstream total_analytic(file6);
total_analytic << "#time   rl    rv    rt" << endl;

//l2 normalised error relative to analytic solution.
fileName file7 = args.rootPath() / args.caseName() / "error_prognostic.dat";
OFstream error_prognostic(file7);
error_prognostic << "#time   rl error    rv error " << endl;

fileName file8 = args.rootPath() / args.caseName() / "error_diagnostic.dat";
OFstream error_diagnostic(file8);
error_diagnostic << "#time   rl error    rv error " << endl;

scalar rl_total;
scalar rv_total;
scalar rt_total;
scalar rl_error;
scalar rv_error;
scalar rt2_total;
const dimensionedScalar Vtot = sum(mesh.V());

#include "writeDiagnostics.H"
