// theta equation
alphaEff = fvc::interpolate(thermophysicalTransport->kappaEff()/Cp);
label nAdvCorr = ucorr < nOuterCorr-1? 1 : 2;
for(label iAdvCorr = 0; iAdvCorr < nAdvCorr; iAdvCorr++)
{
    fvScalarMatrix thetaEqn
    (
        EulerDdt.fvmDdt(rho, thetap)
      + fvm::div(ImEx*aphi, thetap, "div(phi,thetap)")
      + fvc::div((1-ImEx)*aphi, thetap, "div(phi,thetap)")
      + fvc::average(offCentreGW)*rho*w*dThetaadz
      - fvm::laplacian
        (
            offCentre*alphaEff,
            thetap
        )
      - thetapRHS
    );

    thetaEqn.solve();
}

/*solve
(
    EulerDdt.fvmDdt(rho, thetap)
  - fv::FCTadvectDiffuse
    (
        thetap, rho, 
        ImEx*aphi, 
        (1-ImEx)*aphi,
        offCentre*alphaEff,
        -fvc::average(offCentreGW)*rho*w*dThetaadz + thetapRHS,
        nAdvCorr
    )
);

Info << "After full solve thetap goes from " << min(thetap).value()
     << " to " << max(thetap).value() << endl;
*/

thetapf = fvc::interpolate(thetap);
thetaf = thetapf + thetaaf;
theta = thetaa + thetap;
