thetaRHS = -fvc::div((1-offCentre)*phi, theta, "div(phi,theta)")
          + fvc::laplacian
            (
                (1-offCentre)*alphaEff,
                theta
            );

Urhs = - fvc::div((1-offCentre)*phi, U, "div(phi,U)")
       + rho*fvc::weightedReconstruct
        (
            (1-offCentrePG)*
            (gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()),
            0.25
        )
     + fvc::laplacian((1-offCentre)*rhof*fvc::interpolate(turbulence->nuEff()), U);

rhoRHS = -fvc::div((1-offCentre)*phi);

