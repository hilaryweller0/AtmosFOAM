// theta equation

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

if (advectionType == advType::Implicit)
{
    thetaEqn += fvScalarMatrix(fvm::div(aphi, thetap, "div(phi,thetap)"));
}
else
{
    thetaEqn += fvc::div(aphi, thetap, "div(phi,thetap)");
}
thetaEqn.solve();

/*alphaEff = fvc::interpolate(thermophysicalTransport->alphaEff());
thetapSource = -rho*w*(fvc::grad(thetaaf) & kdir);

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
*/
thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
