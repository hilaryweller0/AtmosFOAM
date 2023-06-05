// Implicit solution of the Momentum equation for the turbulence

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevTau(U)
  + 2*rho*(Omega^(U-Ug))
  + rho*muSponge*(spongeT & U)
);
if (impU)
{
    if (!hydrostatic)
    {
        solve
        (
            UEqn == rho
            *(
                fvc::reconstruct
                (
                    gSf
                  - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
                )
              //- muSponge*(spongeT & U)
            )
        );
    }
    else
    {
        volVectorField gradP
        (
            "gradP",
            rho*fvc::reconstruct(Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf())
        );
        gradP.replace(2, scalar(0));
        solve(UEqn == -gradP);
        //-rho*muSponge*(spongeT & U));
    }
}
