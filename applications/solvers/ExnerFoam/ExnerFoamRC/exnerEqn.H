while (pimple.correct())
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
//      + rhorAUf*fvc::ddtCorr(rho, U, phi)
    );

    while (pimple.correctNonOrthogonal())
    {
	fvScalarMatrix ExnerEqn
	(
	    fvm::ddt(Psi, Exner)
	  + fvc::div(phiHbyA)
	  - fvm::laplacian(rhof*rhorAUf*Cp*thetaf, Exner)
	  + fvc::div(rhof*rhorAUf*gSf)
	);

	ExnerEqn.solve(mesh.solver(Exner.select(pimple.finalInnerIter())));

	if (pimple.finalNonOrthogonalIter())
	{
            phi = phiHbyA - rhof*rhorAUf*Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf()
                + rhof*rhorAUf*gSf;
            U = HbyA + fvc::reconstruct
            (
                (- Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf() + gSf)*rhorAUf
            );
	    Uf = fvc::interpolate(U);
	}
    }
}

