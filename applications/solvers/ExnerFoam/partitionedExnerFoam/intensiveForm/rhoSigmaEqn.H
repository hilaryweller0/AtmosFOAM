// update sigma according to the continuity equation
sumSigma == 0;
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] = 
    (
        sigma[ip]*rho[ip].oldTime()
      - dt*fvc::div(flux[ip], sigma[ip])
    )/rho[ip];
    sumSigma += sigma[ip];
}

Info << "1-sumSigma goes from " << min(1-sumSigma).value() << " to "
     << max(1-sumSigma).value() << endl;

rhoSum = sigma[0]*rho[0];
for(label ip = 1; ip < nParts; ip++)
{
    rhoSum += sigma[ip]*rho[ip];
}

// Calculate mass tranfer terms
for(label ip = 0; ip < nParts; ip++) { for(label jp=ip+1; jp < nParts; jp++)
{
    massTransfer[ip][jp] = 0.5*Ksigma*
    (
        rho[jp]*fvc::laplacian(sigma[jp])
      - rho[ip]*fvc::laplacian(sigma[ip])
    );
    massTransfer[jp][ip] = -massTransfer[ip][jp];
}}

Info << "Before diffusion sigma[1] goes from " << min(sigma[1]).value()
     << " to " << max(sigma[1]).value() << endl;

// Now update with diffusion
for(label ip = 0; ip < nParts; ip++) { for(label jp=0; jp < nParts; jp++)
{
    if (ip != jp)
    {
        sigma[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp])/rho[ip];
    }
}}

Info << "After diffusion sigma[1] goes from " << min(sigma[1]).value()
     << " to " << max(sigma[1]).value() << endl;

for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip]);
}


