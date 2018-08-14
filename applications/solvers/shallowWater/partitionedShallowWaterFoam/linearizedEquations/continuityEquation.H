// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    sigmah[ip] = sigmah[ip].oldTime() - dt*
    (
        fvc::div(VOLFLUX[ip], sigmah[ip])
      + fvc::div(volFlux[ip], SIGMAH[ip])
      //  fvc::reconstruct(VOLFLUX[ip] * fvc::snGrad(sigmah[ip]))
      //+ SIGMAH[ip] * fvc::div(volFlux[ip])
    );
}
//Update total height.
h = sigmah[0] + sigmah[1] + SIGMAH[0] + SIGMAH[1];


