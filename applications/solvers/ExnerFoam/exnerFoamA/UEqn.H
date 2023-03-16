// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevTau(U)
  + 2*rho*(Omega^(U-Ug))
);
/*solve
(
    UEqn == rho*
    (
        fvc::reconstruct
        (
            gSf
          - Cp*fvc::interpolate(theta)*fvc::snGrad(Exner)*mesh.magSf()
        )
      - muSponge*(spongeT & U)
    )
);
*/

