// update rhoSigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRho.ddt()[ip] = -fvc::div(volFlux[ip], sigmaRho[ip], "sigmaRho");

    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*
    (
        (1-offCentre)*sigmaRho.ddt()[ip].oldTime()
      + offCentre*sigmaRho.ddt()[ip]
    );
    
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRho.updateSum();
sigmaRhof.updateSum();

/*// Calculate mass tranfer terms
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
       massTransfer[ip][jp] = Ksigma*
        (
            rho[jp]*fvc::laplacian(sigma[jp])
          - rho[ip]*fvc::laplacian(sigma[ip])
        );
    massTransfer[jp][ip] = -massTransfer[ip][jp];
    }
}

Info << "Before diffusion sigma[0] goes from " << min(sigmaRho[1]/rho[1]).value()
     << " to " << max(sigmaRho[1]/rho[1]).value() << endl;

// Now update sigmaRho with diffusion
for(label ip = 0; ip < nParts; ip++) { for(label jp=0; jp < nParts; jp++)
{
    if (ip != jp)
    {
        sigmaRho[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
    }
}}

Info << "After diffusion sigma[0] goes from " << min(sigmaRho[1]/rho[1]).value()
     << " to " << max(sigmaRho[1]/rho[1]).value() << endl;
*/
