surfaceScalarField rhorAUf
(
    "rhorAUf",
    rhof*
    (
        fvc::interpolate(rho/UEqn.A())/*(1-mag(gUnitNormal))
      + fvc::interpolate(rho/wEqn.A())*mag(gUnitNormal)*/
    )
);

surfaceScalarField G("G", 1 + adt*mag(gUnitNormal)*fvc::interpolate(muSponge));

/*if (SIgravityWaves)
{
}
else
{*/
    G == 1/G;
//}

for (int corr=0; corr<nCorr; corr++)
{
    // Velocity without pressure gradient or buoyancy
    U = UEqn.H()/UEqn.A();
    //w = wEqn.H()/wEqn.A();
    
    // Flux without pressure gradient
    //phi = fvc::flux(rho*(U - w*ghat)) - rhorAUf*gSf*fvc::interpolate(thetap/thetaa);
    phi = G*(fvc::flux(rho*U) - rhorAUf*gSf*fvc::interpolate(thetap/thetaa));

    // Exner from eqn of state
    Exnerp.oldTime() = pow(rho*R*(thetaa+thetap)/pRef, RbyCv) - Exnera;
    Exnerp.oldTime().correctBoundaryConditions();
    volScalarField Psi = rho/(RbyCv*(Exnera+Exnerp.oldTime()));

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix ExnerEqn
        (
            EulerDdt.fvmDdt(Psi, Exnerp)
          + EulerDdt.fvcDdt(rho)
          + (1-ocAlpha)*fvc::div(phi.oldTime())
          + ocAlpha*fvc::div(phi)
          - ocAlpha*fvm::laplacian(G*rhorAUf*Cp*thetaf, Exnerp)
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

            U -= rho*
            (
                g*thetap/thetaa
              + Cp*(thetaa + thetap)*fvc::grad(Exnerp)
            )/UEqn.A();

            U = inv(tensor::I + adt*muSponge*(ghat*ghat)) & U;
                    
            /*volVectorField Urhs = rho*
            (
                g*thetap/thetaa
              - Cp*(thetaa + thetap)*fvc::grad(Exnerp)
            );
            volScalarField wrhs = - (Urhs & ghat);
            Urhs += wrhs*ghat;
            
            U += Urhs/UEqn.A();
            w += wrhs/wEqn.A();*/
        }
    }
}

