volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(sqr(rho)*rAU));

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
              + fvc::div(phi)
              - fvm::laplacian(rhorAUf*Cp*thetaf, Exner)
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
                phi -= rhorAUf*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf();
                U += fvc::reconstruct
                (
                    rhorAUf/rhof*(gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf())
                );
            }
        }
        // Update the pressure and temperature based on the new Exner
        thermo.p() = pRef*pow(Exner, 1/kappa);
        thermo.T() = theta*Exner;
        thermo.he() == thermo.he(thermo.p(),thermo.T());
        thermo.correct();
    }
}

