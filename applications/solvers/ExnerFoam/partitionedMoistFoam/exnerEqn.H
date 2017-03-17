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
        
        surfaceScalarField gradPcoeff2("gradPcoeff2", offCentre*dt*gradPcoeff);
        gradPcoeff2.write();
        volScalarField divU("divU", fvc::div(atmosParts.flux()));
        divU.write();
        atmosParts.flux().write();
        Psi.write();
        Exner.write();
        volVectorField uTmp("uTmp", fvc::reconstruct(atmosParts.flux()));
        uTmp.write();
        FatalErrorIn("") << exit(FatalError);

        // Partition back-substitutions
        for(label ip = 0; ip < atmosParts.size(); ip++)
        {
            partition& parti = atmosParts[ip];
            parti.flux() += -offCentre*dt*parti.gradPcoeff()
                               *fvc::snGrad(Exner)*mesh.magSf()
                        + offCentre*dt*fvc::interpolate(parti.sigmaRho())*gSf;
            parti.u() = fvc::reconstruct
            (
                parti.flux()/fvc::interpolate(parti.sigmaRho())
            );
            parti.Uf() = fvc::interpolate(parti.u());
            parti.dFluxdt() += fvc::interpolate(parti.sigmaRho())*gSf
                          - parti.gradPcoeff()*fvc::snGrad(Exner)*mesh.magSf();
        }
    }
    atmosParts.dRhodt() = -fvc::div(atmosParts.flux());
}

