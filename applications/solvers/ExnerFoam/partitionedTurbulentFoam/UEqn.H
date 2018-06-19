// Implicit solution of the Momentum equations for the turbulence without
// transfer terms

for(label ip = 0; ip < nParts; ip++)
{
    fvVectorMatrix UEqn
    (
        fvm::ddt(u[ip])
      + fvm::div(volFlux[ip], u[ip], "u")
      - fvm::Sp(fvc::div(volFlux[ip]), u[ip])
      + (1/sigmaRho[ip])*turbulence[ip].divDevRhoReff(u[ip])
      + 2*(Omega^u[ip])
    );
    
    solve
    (
        UEqn ==
        fvc::reconstruct
        (
            gSf - Cp*fvc::interpolate(theta[ip])*fvc::snGrad(Exner)*mesh.magSf()
        )
    );

    Uf[ip] = fvc::interpolate(u[ip]);
}

