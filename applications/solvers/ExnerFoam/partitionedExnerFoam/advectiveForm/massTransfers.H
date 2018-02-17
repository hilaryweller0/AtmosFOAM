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


// Mass transfer based on vertical velocity and horizontal divergence
if (wTransfer)
{
    if (nParts != 2)
    {
        FatalErrorIn("rhoSigmaEqn.H")
            << " mass transfer terms based on w formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }
    
    massTransfer[0][1] -= max(sign(u[0].component(vector::Z)), dimlessZero)*
    min
    (
        fvc::div(flux[0] - sigmaRhof[0]*((Uf[0] & ghat)*ghatSf)),
        massTransferZero
    );
    
    massTransfer[1][0] -= min(sign(u[1].component(vector::Z)), dimlessZero)*
    max
    (
        fvc::div(flux[1] - sigmaRhof[1]*((Uf[1] & ghat)*ghatSf)),
        massTransferZero
    );
     volScalarField transfer
     (
        "transfer",
        massTransfer[0][1] - massTransfer[1][0]
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

