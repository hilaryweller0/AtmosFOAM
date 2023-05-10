// flux form rho theta equation
solve
(
    fvm::ddt(rho, theta)
  + (1-offCentre)*divUtheta.oldTime()
  + offCentre*fvc::div(U,theta)
  //- fvm::laplacian(diffuseTheta, theta)
);

