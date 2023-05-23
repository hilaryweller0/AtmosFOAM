surfaceScalarField rhorAUf
(
    "rhorAUf",
    rhof*
    (
        fvc::interpolate(rho/uEqn.A())*(1-mag(gUnitNormal))
      + fvc::interpolate(rho/wEqn.A())*mag(gUnitNormal)
    )
);

// Exner from eqn of state
Exnerp.oldTime() = pow(rho*R*theta/pRef, RbyCv) - Exnera;
Exnerp.oldTime().correctBoundaryConditions();
volScalarField Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) - rhorAUf*gSf*thetapf/thetaaf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi/ocAlpha, Exnerp)
          + EulerDdt.fvcDdt(rho)/ocAlpha
          + (1-ocAlpha)/ocAlpha*fvc::div(phi.oldTime())
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf*Cp*thetaf, Exnerp)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr&&corr==nCorr-1)
        {
            ExnerEqn.solve(Exnerp.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exnerp.name());
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += ExnerEqn.flux();

            volVectorField Urhs = rho*fvc::reconstruct
            (
                ExnerEqn.flux()/rhorAUf - gSf*thetapf/thetaaf
            );

            u += (Urhs & idir)/uEqn.A();
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + w*kdir;

