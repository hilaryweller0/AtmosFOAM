// theta equation
alphaEff = fvc::interpolate(thermophysicalTransport->alphaEff());
thetapSource = -rho*w*(fvc::grad(thetaaf) & kdir);

//volScalarField thetap2("thetap2", thetap);
//thetap2.oldTime() == thetap.oldTime();

fv::FCTadvectionDiffusion
(
    thetap,
    rho,
    offCentre,
    phi, 
    alphaEff,
    thetapSource,
    advectionType == advType::Implicit,
    ucorr == nOuterCorr-1
);
/*
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + fvc::average(offCentre)*rho*w*(fvc::grad(thetaaf) & kdir)
  - fvm::laplacian
    (
        offCentre*fvc::interpolate(thermophysicalTransport->alphaEff()),
        thetap
    )
  - thetapRHS
);

if (advectionType == advType::Implicit || advectionType == advType::Diagonal)
{
    thetaEqn += fvScalarMatrix(fvm::div(aphi, thetap, "div(phi,thetap)"));
}
else
{
    thetaEqn += fvc::div(aphi, thetap, "div(phi,thetap)");
}

if (advectionType == advType::Implicit)
{
    thetaEqn.solve();
}
else
{
    thetap = thetaEqn.H()/thetaEqn.A();
}

thetap2.write();
Info << "thetapDiff goes from " << min(thetap - thetap2).value()
     << " to "  << max(thetap - thetap2).value() << endl;
*/
thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
