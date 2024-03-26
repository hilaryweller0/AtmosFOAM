surfaceScalarField rhorAUf
(
    "rhorAUf",
    rhof*
    (
        fvc::interpolate(rho/uEqn.A())*(1-mag(gUnitNormal))
      + fvc::interpolate(rho/vEqn.A())*(1-mag(gUnitNormal))
      + fvc::interpolate(rho/wEqn.A())*mag(gUnitNormal)
    )
);

// Exner from eqn of state
Exnerp.oldTime() = pow(rho*R*theta/pRef, RbyCv) - Exnera;
Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    v = vEqn.H()/vEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + v*jdir + w*kdir;

    // Flux without pressure gradient
    phi = fvc::flux(rho*U) - offCentrePG*rhorAUf*gSf*thetapf/thetaaf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi, Exnerp)
          + EulerDdt.fvcDdt(rho)
          - rhoRHS.oldTime()
          + fvc::div(offCentre*phi)
          - fvm::laplacian(offCentre*offCentrePG*rhorAUf*Cp*thetaf, Exnerp)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr==nCorr-1)
        {
            ExnerEqn.solve(Exnerp.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exnerp.name());
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += ExnerEqn.flux()/offCentre;

            Urhs = rho*fvc::weightedReconstruct
            (
                offCentrePG*
                (
                    ExnerEqn.flux()/(offCentre*offCentrePG*rhorAUf)
                  - gSf*thetapf/thetaaf
                ),
                0.25
            );

            u += (Urhs & idir)/uEqn.A();
            v += (Urhs & jdir)/vEqn.A();
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + v*jdir + w*kdir;
aphi = offCentre*phi;

