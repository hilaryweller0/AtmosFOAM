forAll(rt_analytic,celli)
{
    //Use the elapsed time to calculate the angle through which
    //the moisture field should have rotated and apply (reverse) rotation.
    rotationAngle = 2*runTime.time().value()*angularVelocity.y();
    d0 = mesh.C()[celli] - rotationCentre;
    d[0] = d0.x()*Foam::cos(rotationAngle) - d0.z()*Foam::sin(rotationAngle);
    d[1] = 0;
    d[2] = d0.z()*Foam::cos(rotationAngle) + d0.x()*Foam::sin(rotationAngle);
    d += rotationCentre;
    d -= distributionCentre;
    
    r = mag(vector(d.x()/halfWidth.x(), 0, d.z()/halfWidth.z()));

    //Assign the rotated data point its initial value
    if (r <= 1)
    {
        rt_analytic[celli] = distributionMax*sqr(Foam::cos(M_PI*r/2));
    }
    else
    {
        rt_analytic[celli] = 0;
    }
    
    rv_analytic = min(rt_analytic,rvs);
    rl_analytic = rt_analytic - rv_analytic;
}
