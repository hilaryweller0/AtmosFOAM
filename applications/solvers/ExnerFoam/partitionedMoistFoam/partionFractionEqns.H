// Advect partition fraction equations for omega
for(label ip = 1; ip < atmosParts.size(); ip++)
{
    partition& part = atmosParts[ip];
    
    part.drhoSigmadt() = -fvc::div(part.flux(), part.sigma());
    part.sigma() = 
    (
        part.sigma().oldTime()*part.rho().oldTime() + dt*
        (
            (1-offCentre)*part.drhoSigmadt().oldTime()
          + offCentre*part.drhoSigmadt()
        )
    )/part.rho();
}


