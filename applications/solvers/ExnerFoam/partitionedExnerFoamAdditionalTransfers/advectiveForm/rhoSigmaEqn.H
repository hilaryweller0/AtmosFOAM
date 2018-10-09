// update rhoSigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    Info << "Min sigmaRho before " << ip << ": " << min(sigmaRho[ip]).value() << endl;
    // Calculate the rate of change using old sigmaRho for bounded advection
    flux[ip] = linearInterpolate(rho[ip])*volFlux[ip];
    
    sigmaRho.ddt()[ip] = -fvc::div
    (
        (1-offCentre)*flux[ip].oldTime() + offCentre*flux[ip],
        sigma[ip].oldTime(),
        "sigmaRho"
    );

    // Update sigmaRho based on the rate of change
    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*sigmaRho.ddt()[ip];
    sigmaRho[ip].correctBoundaryConditions();
    
    Info << "Min sigmaRho after " << ip << ": " << min(sigmaRho[ip]).value() << endl;
    
    // Update the sums and face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
flux.updateSum();
sigmaRho.updateSum();
sigmaRhof.updateSum();

