thetapRHS = -fvc::div((1-offCentre)*phi, thetap, "div(phi,thetap)")
          - (1-fvc::average(offCentrePG))*rho*w*dThetaadz
          + fvc::laplacian
            (
                (1-offCentre)*fvc::interpolate(thermophysicalTransport->alphaEff()),
                thetap
            );

Urhs = - fvc::div((1-offCentre)*phi, U, "div(phi,U)")
       - rho*fvc::weightedReconstruct
        (
            (1-offCentrePG)*
            (gSf*thetapf/thetaaf + Cp*thetaaf*fvc::snGrad(Exnerp)*mesh.magSf()),
            0.25
        )
     + fvc::laplacian((1-offCentre)*rhof*fvc::interpolate(turbulence->nuEff()), U);

rhoRHS = -fvc::div((1-offCentre)*phi);

