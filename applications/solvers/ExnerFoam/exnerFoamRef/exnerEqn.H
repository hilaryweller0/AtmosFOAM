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
if (!Boussinesq)
{
    Exnerp.oldTime() = pow(rho*R*theta/pRef, RbyCv) - Exnera;
    Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));
    //Psi.oldTime() = Psi;
}

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) - offCentrePG*rhorAUf*gSf*thetapf/thetaaf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            fvc::div(offCentre*phi)
          - fvm::laplacian(offCentre*offCentrePG*rhorAUf*Cp*thetaf, Exnerp)
        );
        if (!Boussinesq)
        {
            ExnerEqn += fvScalarMatrix
            (
                EulerDdt.fvmDdt(Psi, Exnerp)
              + EulerDdt.fvcDdt(rho)
              - rhoRHS.oldTime()
            );
        }
        else
        {
            ExnerEqn.setReference(pRefCell, Exnerp[pRefCell]);
        }

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
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + w*kdir;
aphi = offCentre*phi;

