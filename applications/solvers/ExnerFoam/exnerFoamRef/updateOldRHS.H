thetapRHS = -fvc::div((1-offCentre)*phi, thetap, "div(phi,thetap)")
          - (1-fvc::average(offCentrePG))*rho*w*dThetaadz
          + fvc::laplacian
            (
                (1-offCentre)*alphaEff,
                thetap
            );

Urhs = - fvc::div((1-offCentre)*phi, u, "div(phi,U)")*idir
       - fvc::div((1-offCentre)*phi, v, "div(phi,U)")*jdir
       - fvc::div((1-offCentre)*phi, w, "div(phi,U)")*kdir
       - rho*fvc::weightedReconstruct
        (
            (1-offCentrePG)*
            (gSf*thetapf/thetaaf + Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()),
            0.25
        )
     + fvc::laplacian((1-offCentre)*rhof*fvc::interpolate(turbulence->nuEff()), U);

rhoRHS = -fvc::div((1-offCentre)*phi);

