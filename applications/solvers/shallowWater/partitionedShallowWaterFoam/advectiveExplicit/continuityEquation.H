// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    h[ip] = h[ip].oldTime() - dt*
    (
        fvc::div(volFlux[ip], h[ip])
    );
    
    volScalarField total_sink   = 0*h[ip];
    volScalarField total_source = 0*h[ip];
    if (useMassTransfer) { total_sink += dt*sourceMag*sink[ip]*h[ip].oldTime(); }
    for (label jp = 0; jp < nParts; jp++)
    {
        if (jp != ip)
        {
            if (useDiffusion)    { total_sink += max(dt*diffusionCoeff*fvc::laplacian(h[jp].oldTime() - h[ip].oldTime()),0*h[ip]); }
            
            total_source = 0*h[ip];
            if (useMassTransfer) { total_source += dt*sourceMag*sink[jp]*h[jp].oldTime(); }
            if (useDiffusion)    { total_source += max(dt*diffusionCoeff*fvc::laplacian(h[ip].oldTime() - h[jp].oldTime()),0*h[ip]); }
            
            total_source = min(total_source, h[jp].oldTime());
            
            h[ip] += total_source;
        }
    }
    total_sink = min(total_sink, h[ip].oldTime());
    h[ip] += -total_sink;
    
    //Calculate total height.
    if (ip == 0) hSum  = h[ip];
    else         hSum += h[ip];
}
