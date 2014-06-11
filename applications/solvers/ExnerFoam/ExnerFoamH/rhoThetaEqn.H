// update density according to the continuity equation
solve
(
    fvm::ddt(rho) + (1-offCentre)*divU + offCentre*fvc::div(U)
);
rhof = fvc::interpolate(rho);

//vperp = rhof*(Uf - (Uf & H.delta())*H.delta()/sqr(H.magd()));
//thetaf = thetaf.oldTime()
//       - offCentre*dt*(vperp & fvc::interpolate(fvc::grad(theta)))/rhof
//       - offCentre*dt*V*fvc::snGrad(theta)/rhof/H.magd();

// theta equation
fvScalarMatrix thetaEqn
(
    fvm::ddt(rho, theta)
  + (1-offCentre)*divUtheta.oldTime()
  + offCentre*fvc::div(U, theta)
  //- fvm::laplacian(diffuseTheta, theta)
);
thetaEqn.solve();

