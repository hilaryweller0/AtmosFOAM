// Updating diagnostic variables in each partition
// and check that sum(sigma) = 1
volScalarField sumSigma
(
    IOobject("sumSigma", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("sumSigma", dimless, scalar(0))
);

//// Global Exner from eqn of state
//volScalarField ExnerS("ExnerS", pow(R*rhoTheta/pRef, kappa/(1-kappa)));

for(label ip = 0; ip < nParts; ip++)
{
    // Update the density from the eqn of state
    //rho[ip] = pRef*pow(ExnerS, (1-kappa)/kappa)/(R*theta[ip]);
    rho[ip] = rhoTheta/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
    
    // Check that sum(sigma) = 1
    sumSigma += sigma[ip];
    
    // Update partition velocities
    u[ip] = fvc::reconstruct
    (
        flux[ip]/fvc::interpolate(sigmaRho[ip])
    );
    Uf[ip] = fvc::interpolate(u[ip]);
    dFluxdt[ip] += fvc::interpolate(sigmaRho[ip])*gSf
                - gradPcoeff[ip]*fvc::snGrad(Exner)*mesh.magSf();
}

Info << "1-sumSigma goes from " << 1-max(sumSigma).value() << " to "
     << 1-min(sumSigma).value() << endl;

