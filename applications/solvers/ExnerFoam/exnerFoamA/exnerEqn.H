volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", rhof*fvc::interpolate(rho*rAU));
thetaf = fvc::interpolate(theta);

if (SIgravityWaves)
{
    // Apply this only where grad Exner * grad theta <0
    surfaceScalarField gradPT = fvc::snGrad(Exner)*fvc::snGrad(theta);
    applySIg = max(-sign(gradPT), scalar(0))*mag(gUnitNormal);
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
    G = 1/(1 + rhorAUf*muSpongef/rhof - applySIg*sqr(rhof)*Cp*gradPT
                                /fvc::interpolate(UEqn.A()*thetaEqn.A()));
}
else
{
    G = 1/(1 + rhorAUf*muSpongef/rhof);
}

{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));

    for (int corr=0; corr<nCorr; corr++)
    {
        // Velocity without pressure gradient or buoyancy
        U = UEqn.H()/UEqn.A();
        
        // Flux without pressure gradient
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
                fvm::ddt(Psi, Exner)
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
                //phi = G*(phi - rhorAUf*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf());
                phi = G*phi + ExnerEqn.flux();
                /*// Final explicit solve for mass conservation
                if(ucorr == nOuterCorr-1 && corr == nCorr-1)
                {
                    solve
                    (
                        fvm::ddt(Psi, Exner) + fvc::div(phi)
                    );
                }*/
                if (SIgravityWaves)
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
                    U = inv(tensor::I + rAU*rho*muSponge*spongeT) & U;
                }
            }
        }
    }
}

