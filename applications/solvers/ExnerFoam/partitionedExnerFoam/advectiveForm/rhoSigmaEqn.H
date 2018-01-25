// Calculate mass tranfer terms
// Diffusion between partitions
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

// Tranfer from 0 to 1 when w>0 and from 1 to 0 when w<0
/*if (wTransfer)
{
    volScalarField w0 = u[0].component(vector::Z);
    volScalarField w1 = u[1].component(vector::Z);
    forAll(u, cellI)
    {
        if (w[cellI] < 0)
        {
            massTransfer[1][0][cellI] = sigmaRho[1][cellI]*0.5/runTime.deltaT();
            massTransfer[0][1][cellI] = 0;
        }
        else
        {
            massTransfer[0][1][cellI] = sigmaRho[0][cellI]*0.5/runTime.deltaT();
            massTransfer[1][0][cellI] = 0;
        }
}*/

// update rhoSigma according to the continuity equation
for(label ip = 0; ip < nParts; ip++)
{
    // Calculate the rate of change
    sigmaRho.ddt()[ip] = -fvc::div(volFlux[ip], sigmaRho[ip], "sigmaRho");

    // Add mass transer terms
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigmaRho.ddt()[ip] += massTransfer[jp][ip] - massTransfer[ip][jp];
        }
    }

    // Update sigmaRho based on the rate of change
    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*
    (
        (1-offCentre)*sigmaRho.ddt()[ip].oldTime()
      + offCentre*sigmaRho.ddt()[ip]
    );
    
    // Update the surface (face) values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
// Update the sums
sigmaRho.updateSum();
sigmaRhof.updateSum();

