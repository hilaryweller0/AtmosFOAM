volScalarField unity
(
    IOobject("unity", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless, scalar(1)),
    "zeroGradient"
);

/*fvScalarMatrix unityEqn
(
    EulerDdt.fvmDdt(rho, unity)
  + fvc::div(phiMid, unity, "div(phi,thetap)")
);
unityEqn.solve("thetaEqn");*/

unity =
(
    unity.oldTime()*rho.oldTime()
  - runTime.deltaT()*fvc::div(phiMid, unity, "test")
)/rho - 1;

unity.correctBoundaryConditions();

Info << "unity-1 goes from " << min(unity).value() << " to " 
     << max(unity).value() << endl;

