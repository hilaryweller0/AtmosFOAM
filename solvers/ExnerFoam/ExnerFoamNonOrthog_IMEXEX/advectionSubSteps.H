// Sub step advection using and RK3 scheme
{
    dimensionedScalar dtsub = dt/nAdvectionSubSteps;

    Info << "Advection sub steps" << flush;
    for(int subStep = 0; subStep < nAdvectionSubSteps; subStep++)
    {
        Info << " " << subStep << flush;
        u = H.reconstructd(V.oldTime()/rhof);
        divUu1 = fvc::interpolate(fvc::div(U,u)) & H.delta();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu1 -= diffusionCoeff*rhof*
                      (H.delta() & fvc::interpolate(fvc::laplacian(u)));
        }
        V.oldTime() -= dtsub*divUu1;
        
        u = H.reconstructd(V.oldTime()/rhof);
        divUu2 = fvc::interpolate(fvc::div(U,u)) & H.delta();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu2 -= diffusionCoeff*rhof*
                      (H.delta() & fvc::interpolate(fvc::laplacian(u)));
        }
        V.oldTime() += dtsub*(0.75*divUu1 - 0.25*divUu2);
        
        u = H.reconstructd(V.oldTime()/rhof);
        divUu3 = fvc::interpolate(fvc::div(U,u)) & H.delta();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu3 -= diffusionCoeff*rhof*
                      (H.delta() & fvc::interpolate(fvc::laplacian(u)));
        }
        V.oldTime() += dtsub*(1/12.*(divUu1 + divUu2) - 2/3.*divUu3);
    }
    
    // update velocity variables
    u = H.reconstructd(V.oldTime()/rhof);
    Uf = fvc::interpolate(u);
    //vperp = rhof*(Uf - (Uf & H.delta())*H.delta()/sqr(H.magd()));

    Info << endl;
}
