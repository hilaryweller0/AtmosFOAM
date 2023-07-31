// theta equation
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, thetap)
  + rho*w*(fvc::grad(offCentre*thetaaf) & kdir)
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

thetapf = fvc::interpolate(thetap);
if (!Boussinesq)
{
    thetaf = thetapf + thetaaf;
    theta = thetaa + thetap;
}
