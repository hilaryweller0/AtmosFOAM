{
    // Global Eqn of state and from each partition
    atmosParts.ExnerFromState(Exner);
    
    for (int corr=0; corr<nCorr; corr++)
    {
        Psi == atmosParts.rho()/Exner;

        // Contributions from each partition
        
        // Diffusivity in each partition
        atmosParts.setGradPcoeff(gradPcoeff);

        // Total flux
        atmosParts.flux()
             == dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0));
        
        // Flux without pressure gradient for each partition
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            surfaceScalarField& fluxi = parti.flux();

            parti.dFluxdt() = -mesh.Sf() & 
            (
                fvc::interpolate(fvc::div(fluxi, parti.u()))
              + (2*fvc::interpolate(parti.sigmaRho())*(Omega^parti.Uf()))
            );
            fluxi = fluxi.oldTime()
                  + (1-offCentre)*dt*parti.dFluxdt().oldTime()
                  + offCentre*dt*parti.dFluxdt();

            atmosParts.flux() += fluxi;
        }

        rhog = fvc::interpolate(atmosParts.rho())*gSf;
        atmosParts.flux() += offCentre*dt*rhog;
    
        // The Helmholtz equation for the global Exner pressure 
        fvScalarMatrix ExnerEqn
        (
            fvm::ddt(Psi, Exner)
          - (1-offCentre)*atmosParts.dRhodt().oldTime()
          + offCentre*fvc::div(atmosParts.flux())
          - fvm::laplacian(sqr(offCentre)*dt*gradPcoeff, Exner)
        );

        if(ucorr == nOuterCorr-1 && corr == nCorr-1)
        {
            ExnerEqn.solve(mesh.solver(Exner.name() + "Final"));
        }
        else
        {
            ExnerEqn.solve(mesh.solver(Exner.name()));
        }

        // Global back-substitutions
        atmosParts.flux() += ExnerEqn.flux()/offCentre;
        atmosParts.Uf() = fvc::interpolate
        (
            fvc::reconstruct
            (
                atmosParts.flux()/fvc::interpolate(atmosParts.rho())
            )
        );
        
        // Partition back-substitutions
        atmosParts[0].flux() = atmosParts.flux();
        for(label ip = 1; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            parti.flux() += -offCentre*dt*parti.gradPcoeff()
                               *fvc::snGrad(Exner)*mesh.magSf()
                        + offCentre*dt*fvc::interpolate(parti.sigmaRho())*gSf;
            atmosParts[0].flux() -= parti.flux();
        }
        /*
        Info << "Checking that partition fluxes sum to the total flux" << endl;
        surfaceScalarField sumFlux = atmosParts.flux();
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            sumFlux -= atmosParts[ip].flux();
        }
        Info << "Minimum discrepancy = " << min(sumFlux).value()
             << " maximum discrepancy = " << max(sumFlux).value() << endl;
        */
        // Reconstruct full velocities and rates of change
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            parti.u() = fvc::reconstruct
            (
                parti.flux()/fvc::interpolate(parti.sigmaRho())
            );
            parti.Uf() = fvc::interpolate(parti.u());
            parti.dFluxdt() += fvc::interpolate(parti.sigmaRho())*gSf
                          - parti.gradPcoeff()*fvc::snGrad(Exner)*mesh.magSf();
        }
    }
}

