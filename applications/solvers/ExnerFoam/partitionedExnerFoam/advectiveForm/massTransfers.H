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


// Mass transfer based on vertical velocity
if (wTransfer)
{
    if (nParts != 2)
    {
        FatalErrorIn("rhoSigmaEqn.H")
            << " mass transfer terms based on w formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }
    
    // Area of the interface between the partitions
    surfaceVectorField gradSigma = fvc::interpolate(fvc::grad(sigma[0]));
    volScalarField interfaceArea
    (
        "interfaceArea",
        0.5*fvc::surfaceSum(mag(mesh.Sf() & unitVector(gradSigma)))
    );
    
    massTransfer[0][1] += max
    (
        rho[0]*u[0].component(vector::Z)*interfaceArea/cellVolume,
        massTransferZero
    );
    
    massTransfer[1][0] += max
    (
        -rho[1]*u[1].component(vector::Z)*interfaceArea/cellVolume,
        massTransferZero
    );
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

