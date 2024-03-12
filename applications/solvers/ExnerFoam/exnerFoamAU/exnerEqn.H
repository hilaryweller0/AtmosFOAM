volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", rhof*fvc::interpolate(rho*rAU));

if (SIgravityWaves)
{
    // Apply this only where grad Exner * grad theta <0
    surfaceScalarField gradPT = fvc::snGrad(Exner)*fvc::snGrad(theta);
    applySIg = max(-sign(gradPT), scalar(0))*mag(gUnitNormal);
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
    G = 1/(1 + rhorAUf*muSpongef/rhof - applySIg*sqr(rhof)*Cp*gradPT
                                /fvc::interpolate(UEqn.A()*thetaEqn.A()));
}
else
{
    G = 1/(1 + rhorAUf*muSpongef/rhof);
}

// Exner from eqn of state
Exner.oldTime() = pow(rho*R*theta/pRef, RbyCv);
Psi = rho/(RbyCv*Exner.oldTime());

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    U = UEqn.H()/UEqn.A();
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) + rhorAUf*gSf;
    
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi/ocAlpha, Exner)
          + EulerDdt.fvcDdt(rho)/ocAlpha
          + (1-ocAlpha)/ocAlpha*fvc::div(phi.oldTime())
          + fvc::div(G*phi)
          - fvm::laplacian(G*rhorAUf*Cp*thetaf, Exner)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
        {
            ExnerEqn.solve(Exner.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exner.name());
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi = G*phi + ExnerEqn.flux();
            if (SIgravityWaves)
            {
                thetaf -= applySIg*phi*fvc::snGrad(theta)
                            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
            }
            U += rAU*rho*fvc::weightedReconstruct
            (
                gSf
              - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf(),
                0.25
            );
        }
    }
}

