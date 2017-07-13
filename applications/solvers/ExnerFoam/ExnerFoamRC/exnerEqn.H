while (piso.correct())
{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);
    
    thetaf = fvc::interpolate(theta);
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, Exner));
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + rhorAUf*fvc::ddtCorr(rho, U, phi)
    );

    while (piso.correctNonOrthogonal())
    {
	fvScalarMatrix ExnerEqn
	(
	    fvm::ddt(Psi, Exner)
	  + fvc::div(phiHbyA)
	  - fvm::laplacian(rhorAUf*Cp*thetaf, Exner)
	  + fvc::div(rhorAUf*gSf)
	);

	ExnerEqn.solve(mesh.solver(Exner.select(piso.finalInnerIter())));

	if (piso.finalNonOrthogonalIter())
	{
            phi = phiHbyA - rhorAUf*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
                + rhorAUf*gSf;
            U = HbyA + rAU*fvc::reconstruct
            (
                - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf() + gSf
            );
	    Uf = fvc::interpolate(U);
	}
    }
}

