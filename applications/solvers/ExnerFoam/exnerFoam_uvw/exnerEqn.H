surfaceScalarField rhorAUf
(
    "rhorAUf",
    rhof*
    (
        fvc::interpolate(rho/UEqn.A())*(1-mag(gUnitNormal))
      + fvc::interpolate(rho/wEqn.A())*mag(gUnitNormal)
    )
);

thetaf = fvc::interpolate(theta);

if (hydrostatic || SIgravityWaves)
{
    // Apply this only where on vertical pointing faces where Gw>0
    surfaceScalarField applySIg = mag(gUnitNormal)*max
    (
        -sign(fvc::snGrad(Exner)*fvc::snGrad(theta)),
        scalar(0)
    );
    
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
}

{
    for (int corr=0; corr<nCorr; corr++)
    {
        // Velocity without pressure gradient or buoyancy
        U = UEqn.H()/UEqn.A();
        w = wEqn.H()/wEqn.A();
        
        // Flux without pressure gradient
        phi = fvc::flux(rho*(U - w*ghat)) + rhorAUf*gSf;

        // Exner from eqn of state
        Exner.oldTime() = pow(rho*R*theta/pRef, RbyCv);
        Exner.oldTime().correctBoundaryConditions();
        volScalarField Psi = rho/(RbyCv*Exner.oldTime());
    
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
                phi = phi + ExnerEqn.flux();

                volVectorField Urhs = rho*fvc::reconstruct
                (
                    gSf
                  - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
                );
                volScalarField wrhs = - (Urhs & ghat);
                Urhs += wrhs*ghat;
                
                U += Urhs/UEqn.A();
                w += wrhs/wEqn.A();
            }
        }
    }
}

