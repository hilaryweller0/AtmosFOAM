surfaceScalarField rhorAUf
(
    "rhorAUf",
    rhof*
    (
        0.5*(mesh.solutionD()[0]+1)
            *fvc::interpolate(rho/uEqn.A())*(1-mag(gUnitNormal))
      + 0.5*(mesh.solutionD()[1]+1)
            *fvc::interpolate(rho/vEqn.A())*(1-mag(gUnitNormal))
      + 0.5*(mesh.solutionD()[2]+1)
           *fvc::interpolate(rho/wEqn.A())*mag(gUnitNormal)
    )
);

// Exner from eqn of state
Exnerp.oldTime() = pow(rho*R*theta/pRef, RbyCv) - Exnera;
Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    if (mesh.solutionD()[0] == 1) u = uEqn.H()/uEqn.A();
    if (mesh.solutionD()[1] == 1) v = vEqn.H()/vEqn.A();
    if (mesh.solutionD()[2] == 1) w = wEqn.H()/wEqn.A();
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

            if (mesh.solutionD()[0] == 1) u += (Urhs & idir)/uEqn.A();
            if (mesh.solutionD()[1] == 1) v += (Urhs & jdir)/vEqn.A();
            if (mesh.solutionD()[2] == 1) w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + v*jdir + w*kdir;
aphi = offCentre*phi;

