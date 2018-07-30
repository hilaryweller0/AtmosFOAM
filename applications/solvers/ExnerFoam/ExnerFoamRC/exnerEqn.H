while (pimple.correct())
{
    Psi == pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);
    
    thetaf = fvc::interpolate(theta);
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", rhof*fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, Exner));
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        rhof*fvc::flux(HbyA)
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
            surfaceScalarField force
                 = gSf - Cp*thetaf*fvc::snGrad(Exner)*mesh.magSf();
            phi = phiHbyA + rhof*rhorAUf*force;
            U = HbyA + fvc::reconstruct(force*rhorAUf);
            Uf = fvc::interpolate(U);
        }
    }
}

