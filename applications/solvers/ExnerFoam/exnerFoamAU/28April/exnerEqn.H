volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", rhof*fvc::interpolate(rho*rAU));
thetaf = fvc::interpolate(theta);

if (hydrostatic || SIgravityWaves)
{
    // Apply this only where grad Exner * grad theta <0
    surfaceScalarField gradPT = fvc::snGrad(Exner)*fvc::snGrad(theta);
    applySIg = max(-sign(gradPT), scalar(0))*mag(gUnitNormal);
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
    if (!hydrostatic)
    {
        G = 1/(1 /*+ rhorAUf*muSpongef/rhof*/ - applySIg*sqr(rhof)*Cp*gradPT
                                /fvc::interpolate(UEqn.A()*thetaEqn.A()));
    }
    else
    {
        G = 1/(1-applySIg /*+ rhorAUf*muSpongef/rhof*/ - applySIg*sqr(rhof)*Cp*gradPT
                                /fvc::interpolate(UEqn.A()*thetaEqn.A()));
    }
}
else
{
    G = 1; ///(1 + rhorAUf*muSpongef/rhof);
}

{
    for (int corr=0; corr<nCorr; corr++)
    {
        // Exner from eqn of state
        Exner.oldTime() = pow(rho*R*theta/pRef, RbyCv);
        Exner.oldTime().correctBoundaryConditions();
        volScalarField Psi = rho/(RbyCv*Exner.oldTime());
    
        // Velocity without pressure gradient or buoyancy
        U = UEqn.H()/UEqn.A();
        
        // Flux without pressure gradient perturbation
        phi = fvc::flux(rho*U) + rhorAUf*gSf;
        // Staggered version
        if (stagger)
        {
            phi += phi.oldTime() - fvc::flux(rho.oldTime()*U.oldTime());
        }
        
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix ExnerEqn
            (
                EulerDdt.fvmDdt(Psi/ocAlpha, Exner)
              + EulerDdt.fvcDdt(rho)/ocAlpha
              + (1-ocAlpha)/ocAlpha*fvc::div(phi.oldTime())
              + fvc::div(G*phi)
              - fvm::laplacian(G*rhorAUf*Cp*thetaf, Exner)
            );

            if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
            {
                ExnerEqn.solve(Exner.name() + "Final");
            }
            else
            {
                ExnerEqn.solve(Exner.name());
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi = G*phi + ExnerEqn.flux();
                if (SIgravityWaves || hydrostatic)
                {
                    thetaf -= applySIg*phi*fvc::snGrad(theta)
                                /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
                }
                if (stagger)
                {
                    // Full staggered version
                    U = fvc::reconstruct(phi)/rho;
                }
                else
                {
                    U += rAU*rho*fvc::reconstruct
                    (
                        gSf
                      - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
                    );
                    //U = inv(tensor::I + rAU*rho*muSponge*spongeT) & U;
                }
            }
        }
    }
}

