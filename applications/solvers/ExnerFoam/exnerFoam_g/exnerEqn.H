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
Exnerg.oldTime() = pow(rho*R*theta/pRef, RbyCv) - gh/(Cp*theta);
Exnerg.oldTime().correctBoundaryConditions();
volScalarField Psi = rho/(RbyCv*(Exnerg.oldTime() + gh/(Cp*theta)));

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    u = uEqn.H()/uEqn.A();
    w = wEqn.H()/wEqn.A();
    U = u*idir + w*kdir;
    
    // Flux without pressure gradient
    phi = fvc::flux(rho*U) + rhorAUf*ghGradThetaByTheta*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi, Exnerg)
          + EulerDdt.fvcDdt(rho)
          + (1-ocAlpha)*fvc::div(phi.oldTime())
          + ocAlpha*fvc::div(phi)
          - ocAlpha*fvm::laplacian(rhorAUf*Cp*thetaf, Exnerg)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr&&corr==nCorr-1)
        {
            ExnerEqn.solve(Exnerg.name() + "Final");
        }
        else
        {
            ExnerEqn.solve(Exnerg.name());
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += ExnerEqn.flux()/ocAlpha;

            volVectorField Urhs = rho*fvc::reconstruct
            (
                (ghGradThetaByTheta - Cp*thetaf*fvc::snGrad(Exnerg))
                *mesh.magSf()
            );
            u += (Urhs & idir)/uEqn.A();
            w += (Urhs & kdir)/wEqn.A();
        }
    }
}
U = u*idir + w*kdir;

