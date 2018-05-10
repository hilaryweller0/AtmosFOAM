if (runTime.time().value() > 100)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        sink[0] *= 0;
        momentumSink[ip] *= 0;
        momentumSource[ip] *= 0;
    }
}

// Advect h in each partition
for(label ip = 0; ip < nParts; ip++)
{
    sigmah_temp[ip] = sigmah[ip].oldTime();
    sigmah[ip] = sigmah_temp[ip] - dt*
    (
        (1-offCentre)*fvc::div(volFlux[ip].oldTime(), sigmah[ip].oldTime())
      + offCentre*fvc::div(volFlux[ip], sigmah[ip])
    );
    
    if (useOperatorSplitting)
    {
        sigmah_temp[ip] = sigmah[ip];
    }
}
//Update total height.
sigmah.updateSum();
h = sigmah.sum();


