thetapRHS = -fvc::div((1-offCentre)*phi, thetap, "div(phi,thetap)")
          - rho*w*(fvc::grad((1-offCentre)*thetaaf) & kdir)
          + fvc::laplacian
            (
                (1-offCentre)*fvc::interpolate(thermophysicalTransport->alphaEff()),
                thetap
            );

if (!hydrostatic)
{
    Urhs = - fvc::div((1-offCentre)*phi, U, "div(phi,U)")
           - rho*fvc::weightedReconstruct
            (
               (1-offCentre)
                *(gSf*thetapf/thetaaf + Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()),
                0.25
            )
         + fvc::laplacian((1-offCentre)*rhof*fvc::interpolate(turbulence->nuEff()), U);
}
else
{
    Urhs = - rho*fvc::weightedReconstruct
            (
                (1-offCentre)
               *(gSf*thetapf/thetaaf + Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()),
                0.25
            );
    Urhs += idir*
    (
         -fvc::div((1-offCentre)*phi, u, "div(phi,U)")
        + fvc::laplacian((1-offCentre)*rhof*fvc::interpolate(turbulence->nuEff()), u)
    );
}

rhoRHS = -fvc::div((1-offCentre)*phi);

