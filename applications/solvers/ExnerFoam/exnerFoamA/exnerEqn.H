volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(sqr(rho)*rAU));

if (SIgravityWaves)
{
    // Apply this only where grad Exner * grad theta <0
    surfaceScalarField gradPT = fvc::snGrad(Exner)*fvc::snGrad(theta);
    surfaceScalarField applySIg = max(-sign(gradPT), scalar(0));
    thetaf += applySIg*phi*fvc::snGrad(theta)
            /(fvc::interpolate(thetaEqn.A())*mesh.magSf());
    G = 1/(1 - applySIg*sqr(rhof)*Cp*gradPT
                /fvc::interpolate(UEqn.A()*thetaEqn.A()));
}

{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));

    for (int corr=0; corr<nCorr; corr++)
    {
        // Velocity without pressure gradient or buoyancy
        U = UEqn.H()/UEqn.A();
        
        // Flux without pressure gradient
        phi = fvc::flux(rho*U) + rhorAUf*gSf;
        
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
                phi = G*(phi - rhorAUf*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf());
                U += rAU*rho*(g - Cp*theta*fvc::grad(Exner));
                if(ucorr == nOuterCorr-1 && corr == nCorr-1)
                {
                    solve
                    (
                        fvm::ddt(Psi, Exner) + fvc::div(phi)
                    );
                }
            }
        }
    }
}

