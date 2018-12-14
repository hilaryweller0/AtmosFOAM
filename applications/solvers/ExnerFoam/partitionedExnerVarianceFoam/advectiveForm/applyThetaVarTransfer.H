// Valid only for 2 fluids
const scalar minSigmaDiv = max(minSigma, VSMALL);
const dimensionedScalar thetaMin("",theta[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaVarMin("",thetaVar[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaOne("",theta[0].dimensions(),scalar(1));
const dimensionedScalar thetaVarOne("",thetaVar[0].dimensions(),scalar(1));
const dimensionedScalar transferOne("",thetaVarTransferTerm[0][0].dimensions(),scalar(1));

partitionedVolScalarField sigmaRhoPreTransfer = sigmaRho;
partitionedVolScalarField thetaPreTransfer = theta;
partitionedVolScalarField thetaVarPreTransfer = thetaVar;

for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " pre Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaVar" << ip << " pre Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
}

massTransfer[0][1] = sigmaRho[0]/thetaVarTimescale * 0.5 *
(
    1 + Foam::erf
    (
        (theta[0] - theta[1]) /
        max
        (
            thetaMin, 
            Foam::sqrt( 2 * (thetaVar[0] + thetaVar[1]) ) 
        )
    )
);

massTransfer[1][0] = sigmaRho[1]/thetaVarTimescale * 0.5 *
(
    1 + Foam::erf
    (
        (theta[1] - theta[0]) /
        max
        (
            thetaMin, 
            Foam::sqrt( 2 * (thetaVar[0] + thetaVar[1]) ) 
        )
    )
);

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

    // Update the face values
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip], "sigmaRho");
}
#include "sigma.H"

//Set transfer terms based on mean values of theta
thetaTransferTerm[0][1] = 1./thetaVarTimescale * 0.5 *
(
    1 - Foam::erf
    (
        (theta[1] - theta[0]) / 
        max
        (
            thetaMin, 
            Foam::sqrt( 2 * (thetaVar[0] + thetaVar[1]) ) 
        )
    )
);
thetaTransferTerm[1][0] = thetaTransferTerm[0][1];

for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            thetaVarTransferTerm[ip][jp] = 1./thetaVarTimescale * Foam::sqrt
            (
                thetaVar[ip] / 
                ( 
                    8*M_PI * max( thetaVarMin, thetaVar[ip] + thetaVar[jp] ) 
                )
            ) * 
            Foam::exp
            (
                - 0.5 * sqr(theta[jp] - theta[ip]) / 
                max(thetaVarMin, thetaVar[jp]) *
                ( 
                    scalar(1) - thetaVar[ip] / 
                    max(thetaVarMin, thetaVar[ip] + thetaVar[jp])  
                )
            );
        }
    }
}

//Theta transfer based on mean theta
theta[0] -= dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRhoPreTransfer[1] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));
theta[1] += dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));

theta[0] += dt*thetaTransferTerm[1][0]*thetaPreTransfer[1] *
    sigmaRhoPreTransfer[1] / max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));
theta[1] -= dt*thetaTransferTerm[1][0]*thetaPreTransfer[1];

//Theta transfer based on theta variance
theta[0] -= dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]);
theta[1] += dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));

theta[0] -= dt*thetaVarTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) *
    sigmaRhoPreTransfer[1] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));
theta[1] += dt*thetaVarTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]);


//Theta variance transfer based on mean theta
thetaVar[0] -= dt*thetaTransferTerm[0][1]*thetaVarPreTransfer[0];
thetaVar[1] += dt*thetaTransferTerm[0][1]*thetaVarPreTransfer[0] * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));

thetaVar[0] += dt*thetaTransferTerm[1][0]*thetaVarPreTransfer[1] *
    sigmaRhoPreTransfer[1] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));
thetaVar[1] -= dt*thetaTransferTerm[1][0]*thetaVarPreTransfer[1];

//Theta variance transfer based on theta variance
thetaVar[0] -= dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) *
    (theta[1] - theta[0]) *
    thetaVar[0]/max(thetaVar[0] + thetaVar[1], thetaVarMin);
thetaVar[1] += dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) * 
    (theta[1] - theta[0]) *
    thetaVar[0]/max(thetaVar[0] + thetaVar[1], thetaVarMin) *
    sigmaRhoPreTransfer[0] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));

thetaVar[0] += dt*thetaTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) *
    (theta[1] - theta[0]) *
    thetaVar[1]/max(thetaVar[0] + thetaVar[1], thetaVarMin) *
    sigmaRhoPreTransfer[1] /  max(sigmaRho.sum(), minSigmaDiv*(rho[0]+rho[1]));
thetaVar[1] -= dt*thetaTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) * 
    (theta[1] - theta[0]) *
    thetaVar[1]/max(thetaVar[0] + thetaVar[1], thetaVarMin);


// Update the sum
sigmaRho.updateSum();
theta.updateSum();

for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " post Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
    Info << "Total theta" << ip << ": " << sum(theta.sum()).value();
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaVar" << ip << " post Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
}
