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

        atmosParts.flux().write();
        surfaceScalarField gradPcoeff2("gradPcoeff2", offCentre*dt*gradPcoeff);
        gradPcoeff2.write();
        FatalErrorIn("") << exit(FatalError);

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

        // Partition back-substitutions
        atmosParts.flux()
             == dimensionedScalar("flux", dimensionSet(1,0,-1,0,0), scalar(0));
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            parti.flux() += offCentre*dt*fvc::interpolate(parti.sigmaRho())*
            (
               - air.Cp()*parti.thetaRho()*fvc::snGrad(Exner)*mesh.magSf()
               + gSf
            );
//            parti.flux() += fvc::interpolate(parti.sigma())*ExnerEqn.flux()/offCentre
//                         + offCentre*dt*fvc::interpolate(parti.sigmaRho())*gSf;
            
            atmosParts.flux() += parti.flux();
        }
        
        // Reconstruct full velocities and rates of change
        atmosParts.Uf() == dimensionedVector("Uf", dimVelocity, vector::zero);
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            parti.u() = fvc::reconstruct
            (
                parti.flux()/fvc::interpolate(parti.sigmaRho())
            );
            parti.Uf() = fvc::interpolate(parti.u());
            parti.dFluxdt() += fvc::interpolate(parti.sigmaRho())*
            (
               - air.Cp()*parti.thetaRho()*fvc::snGrad(Exner)*mesh.magSf()
               + gSf
            );
            atmosParts.Uf() += fvc::interpolate(parti.sigma())*parti.Uf();
        }
    }
}

