// sigmaRho on the face for advection and for the Helmholtz equation
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip]);
}
sigmaRhof.updateSum();
