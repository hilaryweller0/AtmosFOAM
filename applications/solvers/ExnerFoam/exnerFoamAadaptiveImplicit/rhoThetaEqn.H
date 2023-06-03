// update density according to the continuity equation
solve
(
    EulerDdt.fvmDdt(rho)
  + fvc::div(phiOld)
  + fvc::div(phiNew)
);
rhof = fvc::interpolate(rho);

// theta equation
fvScalarMatrix thetaEqn
(
    EulerDdt.fvmDdt(rho, theta)
  + fvm::div(phiNew, theta, "div(phi,theta)")
  + fvc::div(phiOld, theta.oldTime(), "div(phi,theta)")
  //- fvm::laplacian(thermophysicalTransport->alphaEff(), theta)

    // Term to compensate energy for numerical diffusion of theta
  /*- 0.5/Exner*fvc::laplacian
    (
        mag(phi)*fvc::interpolate(theta, "diffusionThetaCorr")
                /(mesh.deltaCoeffs()*mesh.magSf()),
        Exner
    )*/
);

if (implicitT)
{
    thetaEqn.solve();
}
else
{
    theta = thetaEqn.H()/thetaEqn.A();
}
thetaf = fvc::interpolate(theta);
