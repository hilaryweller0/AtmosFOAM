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
Exner.oldTime() = pow(rho*R*theta/pRef, RbyCv);
volScalarField Psi = rho/(RbyCv*Exner.oldTime());

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;

    // Flux without pressure gradient
    phi = fvc::flux(rho*U) + rhorAUf*gSf;
    
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi, Exner)
          + EulerDdt.fvcDdt(rho)
          + (1-ocAlpha)*fvc::div(phi.oldTime())
          + ocAlpha*fvc::div(phi)
          - ocAlpha*fvm::laplacian(rhorAUf*Cp*thetaf, Exner)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
        {
            ExnerEqn.solve(Exner.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exner.name());
        }

        if (nonOrth == nNonOrthCorr && ocAlpha > SMALL)
        {
            phi += ExnerEqn.flux()/ocAlpha;

            volVectorField Urhs = rho*fvc::weightedReconstruct
            (
                gSf + ExnerEqn.flux()/(ocAlpha*rhorAUf), 0.25
            );

            u += (Urhs & idir)/uEqn.A();
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + w*kdir;

