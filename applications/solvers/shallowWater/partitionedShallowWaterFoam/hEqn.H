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
                         - fvc::interpolate(fvc::div(flux[ip], u[ip]))
                       ));
            
            fluxSum += flux[ip];
        }

        // Helmholtz equation
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvc::div(fluxSum)
              - fvm::laplacian(dt*g*h, h)
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
            if (nonOrth == nNonOrthCorr) { for(label ip = 0; ip < nParts; ip++)
            {
                surfaceScalarField sigmaf = fvc::interpolate(sigma[ip]);
                flux[ip] -= dt*sigmaf*g*fvc::interpolate(h)*fvc::snGrad(h)*mesh.magSf();
                u[ip] = fvc::reconstruct(flux[ip])/(sigma[ip]*h);
                Uf[ip] = fvc::interpolate(u[ip]);
            }}
        }
    }
}
