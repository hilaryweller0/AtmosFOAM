// update rhoSigma according to the continuity equation (without transfers)
for(label ip = 0; ip < nParts; ip++)
{
    // Calculate the rate of change using old sigmaRho for bounded advection
    sigmaRho.ddt()[ip] = -fvc::div
    (
        (1-offCentre)*volFlux[ip].oldTime() + offCentre*volFlux[ip],
        sigmaRho[ip].oldTime(),
        "sigmaRho"
    );

    // Update sigmaRho based on the rate of change
    sigmaRho[ip] = sigmaRho[ip].oldTime() + dt*sigmaRho.ddt()[ip];
    sigmaRho[ip].correctBoundaryConditions();
    
/*    Tuple2<scalar, label> minSigmaRho = minWhere(sigmaRho[ip]);
    label cellI = minSigmaRho.second();
    Info << "sigmaRho[" << ip << "] minimum " 
         << minSigmaRho.first() << " at " << cellI
         << " where sigmaRho[" << ip << "].oldTime() is " 
         << sigmaRho[ip].oldTime()[cellI]
         << " and dt*sigmaRho.ddt()[" << ip << "] is "
         << dt.value()*sigmaRho.ddt()[ip][cellI] << endl;
*/
}

// Calculate mass tranfer terms as diffusion between partitions
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransferTmp = 0.5*Ksigma*fvc::laplacian
        (
            sigmaRho[jp] -sigmaRho[ip]
        );
        massTransfer[ip][jp] = max(massTransferTmp, massTransferZero);
        massTransfer[jp][ip] = -min(massTransferTmp, massTransferZero);
    }
}


// Mass transfer based on buoyancy
if (thetaTransfer)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        surfaceScalarField gradTheta = fvc::snGrad(theta[ip], "theta");
        gradTheta -= gradTheta*gUnitNormal;
        horizLaplacianTheta[ip] = fvc::div(gradTheta*mesh.magSf());
    }
    if (nParts != 2)
    {
        FatalErrorIn("rhoSigmaEqn.H")
            << " mass transfer terms only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }
    for(label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (horizLaplacianTheta[0][cellI] < 0)
        {
            massTransfer[0][1][cellI] += sigmaRho[0][cellI]/dt.value();
        }
        if (horizLaplacianTheta[1][cellI] > 0)
        {
            massTransfer[1][0][cellI] += sigmaRho[1][cellI]/dt.value();
        }
    }
}

// Constraint to keep sigmaRho > 0
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp], (sigmaRho[ip]-minSigma*rho[ip])/dt
            );
        }
    }
}

/*for(label ip = 0; ip < nParts; ip++)
{
    Tuple2<scalar, label> minSigmaRho = minWhere(sigmaRho[ip]);
    label cellI = minSigmaRho.second();
    Info << "sigmaRho[" << ip << "] minimum " 
         << minSigmaRho.first() << " at " << cellI
         << " where sigmaRho[" << ip << "].oldTime() is " 
         << sigmaRho[ip].oldTime()[cellI]
         << " and dt*sigmaRho.ddt()[" << ip << "] is "
         << dt.value()*sigmaRho.ddt()[ip][cellI] << endl;
}
*/
Info << "Before mass transfers, sigmaRho[0] goes from "
     << min(sigmaRho[0]).value() << " to " << max(sigmaRho[0]).value() << endl;
Info << "Before mass transfers, sigmaRho[1] goes from "
     << min(sigmaRho[1]).value() << " to " << max(sigmaRho[1]).value() << endl;

// Apply mass transfer terms (operator split) to sigmaRho equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigmaRho[ip] += dt*(massTransfer[jp][ip] - massTransfer[ip][jp]);
        }
    }
}

Info << "After mass transfers, sigmaRho[0] goes from "
     << min(sigmaRho[0]).value() << " to " << max(sigmaRho[0]).value() << endl;
Info << "After mass transfers, sigmaRho[1] goes from "
     << min(sigmaRho[1]).value() << " to " << max(sigmaRho[1]).value() << endl;

// Update the sums and face values
for(label ip = 0; ip < nParts; ip++)
{
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
sigmaRho.updateSum();
sigmaRhof.updateSum();

// Transfer terms for other variables (assuming nParts == 2)
if (nParts != 2)
{
    FatalErrorIn("rhoSigmaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}
transfer[0] = dt*massTransfer[0][1]/max(sigmaRho[1], minSigma*rho[1]);
transfer[1] = dt*massTransfer[1][0]/max(sigmaRho[0], minSigma*rho[0]);
transfer.updateSum();
transferf[0] = fvc::interpolate(transfer[0], "massTransfer");
transferf[1] = fvc::interpolate(transfer[1], "massTransfer");
transferf.updateSum();

