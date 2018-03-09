// Calculate mass tranfer terms as diffusion between partitions
if (Ksigma.value() > SMALL)
{
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
}
else
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            massTransfer[ip][jp] == massTransferZero;
            massTransfer[jp][ip] == massTransferZero;
        }
    }
}

//Info << "Before divTransfer"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

// Mass transfer to avoid changes in sigmaRho due to partition divergence
if (divTransfer)
{
    if (nParts != 2)
    {
        FatalErrorIn("rhoSigmaEqn.H")
            << " mass transfer terms based on divTransfer formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }
    
    for(label ip = 0; ip < partNames.size(); ip++)
    {
        divu[ip] = fvc::div(volFlux[ip]);
    }
    
    volScalarField T = sigmaRho[0]*(divu[0] - divu.updateSum()/rho.sum());
    massTransfer[0][1] -= min(T, massTransferZero);
    massTransfer[1][0] += max(T, massTransferZero);
}

//Info << "After divTransfer"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

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
        fvc::div(sigmaFlux[0] - sigmaRhof[0]*((Uf[0] & ghat)*ghatSf)),
        massTransferZero
    );
    
    massTransfer[1][0] -= min(sign(u[1].component(vector::Z)), dimlessZero)*
    max
    (
        fvc::div(sigmaFlux[1] - sigmaRhof[1]*((Uf[1] & ghat)*ghatSf)),
        massTransferZero
    );
}

//Info << "After wTransfer"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

// Mass transfer based on buoyancy
if (thetaTransfer)
{
    if (nParts != 2)
    {
        FatalErrorIn("rhoSigmaEqn.H")
            << " mass transfer terms only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    massTransfer[0][1] -= sigmaRho[0]*thetaTransferDiffusivity*min
    (
        fvc::laplacian(theta[0])/theta[0],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
    massTransfer[1][0] += sigmaRho[1]*thetaTransferDiffusivity*max
    (
        fvc::laplacian(theta[1])/theta[1],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

//Info << "After thetaTransfer"
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value() << endl;

// Constraint to keep sigmaRho > 0
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max
                (
                    (max(sigmaRho[ip], densityZero)-minSigma*rho[ip])/dt,
                    massTransferZero
                )
            );
        }
    }
}

//Info << "Before mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value()
//     << "\nmassTransfer[0][1] goes from " << min(massTransfer[0][1]).value()
//     << " to " << max(massTransfer[0][1]).value()
//     << "\nmassTransfer[1][0] goes from " << min(massTransfer[1][0]).value()
//     << " to " << max(massTransfer[1][0]).value()
//     << "\nmin(rho[0]) = " << min(rho[0]).value() 
//     << "\nmin(rho[1]) = " << min(rho[1]).value() << endl;

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

//Info << "After mass transfer min(sigmaRho[0]) = " << min(sigmaRho[0]).value()
//     << " and min(sigmaRho[1]) = " << min(sigmaRho[1]).value() << endl;

// Update the sum
sigmaRho.updateSum();

// Transfer terms for other variables (assuming nParts == 2)
if (nParts != 2)
{
    FatalErrorIn("rhoSigmaEqn.H")
        << " transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}
const scalar minSigmaDiv = max(minSigma, VSMALL);
transfer[0] = dt*massTransfer[0][1]/max(sigmaRho[1], minSigmaDiv*rho[1]);
transfer[1] = dt*massTransfer[1][0]/max(sigmaRho[0], minSigmaDiv*rho[0]);
transfer.updateSum();

