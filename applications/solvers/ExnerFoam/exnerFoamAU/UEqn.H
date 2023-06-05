// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevTau(U)
  + 2*rho*(Omega^(U-Ug))
  + rho*muSponge*(spongeT & U)
);
if (implicitU)
{
    solve
    (
        UEqn == rho
        *(
            fvc::weightedReconstruct
            (
                gSf
              - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf(),
                0.25
            )
        )
    );
}
