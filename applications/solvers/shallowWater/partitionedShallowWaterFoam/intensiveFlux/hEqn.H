{
    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of velocity for each partition
        fluxSum
             == dimensionedScalar("flux", dimensionSet(0,4,-1,0,0), scalar(0));
        for(label ip = 0; ip < nParts; ip++)
        {
            // Flux without pressure gradient
            flux[ip] = flux[ip].oldTime()
                     - dt*(mesh.Sf() &
                       (
                        - fvc::interpolate(fvc::div(flux[ip], sigma[ip]*u[ip]))
                       ));
            
            fluxSum += sigmaf[ip]*flux[ip];
        }

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvc::div(fluxSum)
              - fvm::laplacian(dt*g*hf, h)
            );

            if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
            {
                hEqn.solve(mesh.solver(h.name() + "Final"));
            }
            else
            {
                hEqn.solve(mesh.solver(h.name()));
            }

            // Back-substitutions
            if (nonOrth == nNonOrthCorr)
            {
                for(label ip = 0; ip < nParts; ip++)
                {
                    flux[ip] -= dt*g*hf*fvc::snGrad(h)*mesh.magSf();
                    u[ip] = fvc::reconstruct(flux[ip])/h;
                    Uf[ip] = fvc::interpolate(u[ip]);
                }
                fluxSum -= dt*g*hf*fvc::snGrad(h)*mesh.magSf();
                
                surfaceScalarField fluxSumCheck = sigmaf[0]*flux[0];
                for(label ip = 1; ip < nParts; ip++)
                {
                    fluxSumCheck += sigmaf[ip]*flux[ip];
                }
                Info << "fluxSumCheck-fluxSum goes from "
                     << min(fluxSumCheck-fluxSum).value() << " to "
                     << max(fluxSumCheck-fluxSum).value() << endl;
            }
            h *= sumSigma;
        }
    }
}
