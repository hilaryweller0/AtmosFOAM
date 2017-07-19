for (int corr=0; corr<nCorr; corr++)
{
    // Explicit part of velocity for each partition
    for(label ip = 0; ip < nParts; ip++)
    {
        // Flux without pressure gradient
        volFlux[ip] = volFlux[ip].oldTime() - dt*
        (
            0.5*fvc::flux(fvc::div(volFlux[ip], u[ip]))
          + ((twoOmegaf^Uf[ip]) & mesh.Sf())
        );
        
        if (ip == 0) fluxSum = hf[ip]*volFlux[ip];
        else fluxSum += hf[ip]*volFlux[ip];
    }

    // Helmholtz equation
    hH = hSum;
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(hH)
          + fvc::div(fluxSum)
          - fvm::laplacian(dt*g*hfSum, hH)
        );

        if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
        {
            hEqn.solve(mesh.solver(hH.name() + "Final"));
        }
        else
        {
            hEqn.solve(mesh.solver(hH.name()));
        }

        // Back-substitutions
        if (nonOrth == nNonOrthCorr)
        {
            for(label ip = 0; ip < nParts; ip++)
            {
                h[ip] *= hH/hSum;
                volFlux[ip] -= dt*g*fvc::snGrad(hH)*mesh.magSf();
                u[ip] = fvc::reconstruct(volFlux[ip]);
                Uf[ip] = fvc::interpolate(u[ip]);
            }
        }
    }
}

