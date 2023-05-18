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
            EulerDdt.fvmDdt(Psi, Exnerp)
          + EulerDdt.fvcDdt(rho)
          + (1-ocAlpha)*fvc::div(phi.oldTime())
          + ocAlpha*fvc::div(phi)
          - ocAlpha*fvm::laplacian(rhorAUf*Cp*thetaf, Exnerp)
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
            phi += ExnerEqn.flux()/ocAlpha;

            volVectorField Urhs = -rho*
                fvc::reconstruct(Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf());
            u += (Urhs & idir)/uEqn.A();
            Urhs -= rho*fvc::reconstruct(thetapf*gSf/thetaaf);
            w += (Urhs & kdir)/wEqn.A();
            U = u*idir + w*kdir;
        }
    }
}

