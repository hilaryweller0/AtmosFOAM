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
Exner.oldTime().correctBoundaryConditions();
volScalarField Psi = rho/(RbyCv*Exner.oldTime());

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) + rhorAUf*Hodge.ddirToFlux(gd)
        - rhorAUf*Cp*thetaf*Hodge.ddirToFluxCorr(fvc::snGrad(Exner)*Hodge.magd());

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi/ocAlpha, Exner)
          + EulerDdt.fvcDdt(rho)/ocAlpha
          + (1-ocAlpha)/ocAlpha*fvc::div(phi.oldTime())
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf*Cp*thetaf, Exner)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr&&corr==nCorr-1)
        {
            ExnerEqn.solve(Exner.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exner.name());
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += ExnerEqn.flux();

            volVectorField Urhs = rho*Hodge.reconstructd
            (
                gd - Cp*thetaf*fvc::snGrad(Exner)*Hodge.magd()
            );

            u += (Urhs & idir)/uEqn.A();
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + w*kdir;

