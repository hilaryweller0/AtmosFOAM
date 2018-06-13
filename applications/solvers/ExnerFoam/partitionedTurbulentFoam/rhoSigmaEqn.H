// update rhoSigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Calculate the rate of change using old sigmaRho for bounded advection
    flux[ip] = linearInterpolate(rho[ip])*volFlux[ip];
    
    solve
    (
        fvm::ddt(sigmaRho[ip]) + fvc::div(flux[ip])
    );
    
    // Update the sums and face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
flux.updateSum();
sigmaRho.updateSum();
sigmaRhof.updateSum();

