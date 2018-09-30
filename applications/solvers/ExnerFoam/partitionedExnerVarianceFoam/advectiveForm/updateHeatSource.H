if ((heatSource) and (runTime.time().value() <= heatSourceDuration.value()))
{
    /*dimensionedScalar timeRatio = runTime.time()/heatSourceDuration;    
    scalar timeFactor = sqr
    ( 
        Foam::sin( M_PI*timeRatio.value() ) 
    );*/

    dimensionedScalar t1 = (runTime.time() - dt/2)/heatSourceDuration;
    dimensionedScalar t2 = (runTime.time() + dt/2)/heatSourceDuration;    
    scalar timeFactor = 0.5*dt.value() + heatSourceDuration.value() / (4*M_PI) *
    (
        Foam::sin(2*M_PI * t1.value())
      - Foam::sin(2*M_PI * t2.value())
    );

    for(label ip = 0; ip < partNames.size(); ip++)
    {
        theta[ip] += timeFactor * thetaSource[ip] / 
        (
            Cv * Exner * max
            (
                sigmaRho[ip], dimensionedScalar("",dimensionSet(1,-3,0,0,0),scalar(1e-16))
            ) 
        );
    }
}
