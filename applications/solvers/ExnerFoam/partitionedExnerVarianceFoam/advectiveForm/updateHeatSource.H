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
        volScalarField deltaTheta = timeFactor * thetaSource[ip] / (Cv*sigmaRho[ip]*Exner);
        forAll(theta[0], celli)
        {
            if (deltaTheta[celli] < scalar(0))
            {
                Info << "NEGATIVE SOURCE, timefactor: " << timeFactor << endl;
                Info << "NEGATIVE SOURCE, sigmarho: " << sigmaRho[ip][celli] << endl;
                Info << "NEGATIVE SOURCE, exner: " << Exner[celli] << endl;
            }
        }

        theta[ip] += timeFactor * thetaSource[ip] / (Cv*sigmaRho[ip]*Exner);
    }
}
