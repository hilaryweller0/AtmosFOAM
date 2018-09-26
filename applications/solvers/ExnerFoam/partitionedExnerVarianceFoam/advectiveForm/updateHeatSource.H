if ((heatSource) and (runTime.time().value() <= heatSourceDuration.value()))
{
    dimensionedScalar timeRatio = runTime.time()/heatSourceDuration;    
    scalar timeFactor = sqr
    ( 
        Foam::sin( M_PI*timeRatio.value() ) 
    );
    
    for(label ip = 0; ip < partNames.size(); ip++)
    {
        theta[ip] += dt*timeFactor*thetaSource[ip];
    }
}
