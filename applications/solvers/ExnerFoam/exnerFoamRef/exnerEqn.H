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
if (thermoDynamic)
{
    Exnerp.oldTime() = pow(rho*R*theta/pRef, RbyCv) - Exnera;
    Exnerp.oldTime().correctBoundaryConditions();
}
volScalarField Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) - ocAlpha*rhorAUf*gSf*thetapf/thetaaf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            fvc::div(phiOld)
          + fvc::div(offCentre*phi)
          - fvm::laplacian(offCentre*ocAlpha*rhorAUf*Cp*thetaf, Exnerp)
        );
        if (thermoDynamic)
        {
            ExnerEqn += fvScalarMatrix
            (
                EulerDdt.fvmDdt(Psi, Exnerp)
              + EulerDdt.fvcDdt(rho)
            );
        }
        else ExnerEqn.relax();

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr&&corr==nCorr-1)
        {
            ExnerEqn.solve(Exnerp.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exnerp.name());
        }

        if (nonOrth == nNonOrthCorr && ocAlpha > SMALL)
        {
            phi += ExnerEqn.flux()/offCentre;

            Urhs = rho*fvc::weightedReconstruct
            (
                ExnerEqn.flux()/(offCentre*ocAlpha*rhorAUf) - gSf*thetapf/thetaaf, 1
            );

            u += ocAlpha*(Urhs & idir)/uEqn.A();
            w += ocAlpha*(Urhs & kdir)/wEqn.A();
        }
    }
}
Urhs.relax(0.5);
U = u*idir + w*kdir;

