// Implicit solution of the Momentum equations for the turbulence without
// transfer terms

for(label ip = 0; ip < nParts; ip++)
{
    fvVectorMatrix UEqn
    (
        fvm::ddt(sigmaRho[ip], u[ip])
      + fvm::div(sigmaFlux[ip], u[ip], "u")
      + sigma[ip]*turbulence[ip].divDevRhoReff(u[ip])
      + 2*sigmaRho[ip]*(Omega^u[ip])
    );
    
    solve
    (
        UEqn ==
        fvc::reconstruct(sigmaRhof[ip]*
        (
            gSf - Cp*fvc::interpolate(theta[ip])*fvc::snGrad(Exner)*mesh.magSf()
        ))
    );

    Uf[ip] = fvc::interpolate(u[ip]);
}
