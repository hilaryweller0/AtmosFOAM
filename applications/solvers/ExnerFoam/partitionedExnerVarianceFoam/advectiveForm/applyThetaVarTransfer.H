Info << "Theta transfers based on theta variance:" << endl;

// Valid only for 2 fluids
const scalar minSigmaDiv = max(minSigma, VSMALL);
const dimensionedScalar thetaMin("",theta[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaVarMin("",thetaVar[0].dimensions(),scalar(1e-16));
const dimensionedScalar thetaOne("",theta[0].dimensions(),scalar(1));
const dimensionedScalar thetaVarOne("",thetaVar[0].dimensions(),scalar(1));
const dimensionedScalar transferOne("",thetaVarTransferTerm[0][0].dimensions(),scalar(1));

thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);

partitionedVolScalarField sigmaRhoPreTransfer = sigmaRho;
partitionedVolScalarField thetaPreTransfer = theta;
partitionedVolScalarField thetaVarPreTransfer = thetaVar;
partitionedVolScalarField thetaNormalised = theta;
partitionedVolScalarField deltaTheta = theta;

for(label ip = 0; ip < nParts; ip++)
{
    //Info << "theta" << ip << " pre Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    //Info << "thetaVar" << ip << " pre Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
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
        (theta[0] - theta[1]) /
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
sigmaRho.updateSum();
#include "sigma.H"



/*transfer[0] += dt*massTransfer[0][1]/max(sigmaRho[1], minSigmaDiv*rho[1]);
transfer[1] += dt*massTransfer[1][0]/max(sigmaRho[0], minSigmaDiv*rho[0]);
// Include the mass transfer and drag implicitly (for nParts=2)
transferf[0] = fvc::interpolate(transfer[0], "massTransfer");
transferf[1] = fvc::interpolate(transfer[1], "massTransfer");
transferf.updateSum();

volFlux[0] =
(
    (1+transferf[0])*volFlux[0] + transferf[1]*volFlux[1]
)/(1+transferf.sum());
volFlux[1] = (volFlux[1] + transferf[0]*volFlux[0])
            /(1+transferf[0]);

// Update mass flux and velocities based on the flux
for(label ip = 0; ip < nParts; ip++)
{
    sigmaFlux[ip] = sigmaRhof[ip]*volFlux[ip];
    u[ip] = fvc::reconstruct(volFlux[ip]);
    Uf[ip] = fvc::interpolate(u[ip], "Uf");
}
sigmaFlux.updateSum();
volFlux.updateSum();*/

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

//Correction to remove transfer if no initial mass present
thetaTransferTerm[0][1] *= sigmaRhoPreTransfer[0] /  
    max(sigmaRhoPreTransfer[0], minSigmaDiv*rho[0]);
thetaTransferTerm[1][0] *= sigmaRhoPreTransfer[1] /  
    max(sigmaRhoPreTransfer[1], minSigmaDiv*rho[1]);
thetaVarTransferTerm[0][1] *= sigmaRhoPreTransfer[0] /  
    max(sigmaRhoPreTransfer[0], minSigmaDiv*rho[0]);
thetaVarTransferTerm[1][0] *= sigmaRhoPreTransfer[1] /  
    max(sigmaRhoPreTransfer[1], minSigmaDiv*rho[1]);

//Normalise theta for internal energy conservation.
thetaNormalised[0] = thetaPreTransfer[0] * 
    sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaNormalised[1] = thetaPreTransfer[1] * 
    sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);
    
theta[0] = thetaNormalised[0];
theta[1] = thetaNormalised[1];
thetaVar[0] *= sigmaRhoPreTransfer[0] / max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[1] *= sigmaRhoPreTransfer[1] / max(sigmaRho[1], minSigmaDiv*rho[1]);

//Theta transfer based on mean theta
deltaTheta[0] = -dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
deltaTheta[1] =  dt*thetaTransferTerm[0][1]*thetaPreTransfer[0] * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho[1], minSigmaDiv*rho[1]);
    
deltaTheta[0] -= dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
deltaTheta[1] += dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) * 
    sigmaRhoPreTransfer[0] /  max(sigmaRho[1], minSigmaDiv*rho[1]);
    
theta[0] += deltaTheta[0];
theta[1] += deltaTheta[1];

thetaVar[0] -= dt*thetaTransferTerm[0][1]*
(
    thetaVarPreTransfer[0] 
  + sqr(thetaNormalised[0]+deltaTheta[0]-thetaPreTransfer[0])
) * sigmaRhoPreTransfer[0] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[1] += dt*thetaTransferTerm[0][1]*
(
    thetaVarPreTransfer[0] 
  + sqr(thetaNormalised[0]+deltaTheta[0]-thetaPreTransfer[0])
) * sigmaRhoPreTransfer[0] /  max(sigmaRho[1], minSigmaDiv*rho[1]);

//Theta variance transfer based on theta variance
thetaVar[0] -= dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) *
(
    thetaPreTransfer[0] - (thetaNormalised[0]+deltaTheta[0])
  + (thetaPreTransfer[1] - thetaPreTransfer[0]) *
    thetaVarPreTransfer[0]/max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) 
) * sigmaRhoPreTransfer[0] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[1] += dt*thetaVarTransferTerm[0][1]*sqrt(thetaVarPreTransfer[0]) * 
(
    thetaPreTransfer[0] - (thetaNormalised[0]+deltaTheta[0])
  + (thetaPreTransfer[1] - thetaPreTransfer[0]) *
    thetaVarPreTransfer[0]/max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin) 
) * sigmaRhoPreTransfer[0] /  max(sigmaRho[1], minSigmaDiv*rho[1]);

deltaTheta[0] =  dt*thetaTransferTerm[1][0]*thetaPreTransfer[1] *
    sigmaRhoPreTransfer[1] / max(sigmaRho[0], minSigmaDiv*rho[0]);
deltaTheta[1] = -dt*thetaTransferTerm[1][0]*thetaPreTransfer[1] * 
    sigmaRhoPreTransfer[1] /  max(sigmaRho[1], minSigmaDiv*rho[1]);
    
deltaTheta[0] -= dt*thetaVarTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) *
    sigmaRhoPreTransfer[1] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
deltaTheta[1] += dt*thetaVarTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) * 
    sigmaRhoPreTransfer[1] /  max(sigmaRho[1], minSigmaDiv*rho[1]);
    
theta[0] += deltaTheta[0];
theta[1] += deltaTheta[1];

thetaVar[0] += dt*thetaTransferTerm[1][0]*
(
    thetaVarPreTransfer[1] 
  + sqr(thetaNormalised[1]+deltaTheta[1]-thetaPreTransfer[1])
) * sigmaRhoPreTransfer[1] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[1] -= dt*thetaTransferTerm[1][0]*
(
    thetaVarPreTransfer[1] 
  + sqr(thetaNormalised[1]+deltaTheta[1]-thetaPreTransfer[1])
) * sigmaRhoPreTransfer[1] /  max(sigmaRho[1], minSigmaDiv*rho[1]);


thetaVar[0] += dt*thetaTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) *
(
    thetaPreTransfer[1] - (thetaNormalised[1]+deltaTheta[1])
  + (thetaPreTransfer[1] - thetaPreTransfer[0]) *
    thetaVarPreTransfer[1]/max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin)
) * sigmaRhoPreTransfer[1] /  max(sigmaRho[0], minSigmaDiv*rho[0]);
thetaVar[1] -= dt*thetaTransferTerm[1][0]*sqrt(thetaVarPreTransfer[1]) * 
(
    thetaPreTransfer[1] - (thetaNormalised[1]+deltaTheta[1])
  + (thetaPreTransfer[1] - thetaPreTransfer[0]) *
    thetaVarPreTransfer[1]/max(thetaVarPreTransfer[0] + thetaVarPreTransfer[1], thetaVarMin)
) * sigmaRhoPreTransfer[1] /  max(sigmaRho[1], minSigmaDiv*rho[1]);

// Update the sum
theta.updateSum();
//thetaVar[0] = max(thetaVar[0], 0*thetaVar[0]);
//thetaVar[1] = max(thetaVar[1], 0*thetaVar[1]);
thetaVar[0] = mag(thetaVar[0]);
thetaVar[1] = mag(thetaVar[1]);

Info << "thetaTransfer01: [" << min(thetaTransferTerm[0][1]).value() << ", " << max(thetaTransferTerm[0][1]).value() << "]" << endl;
Info << "thetaTransfer10: [" << min(thetaTransferTerm[1][0]).value() << ", " << max(thetaTransferTerm[1][0]).value() << "]" << endl;

for(label ip = 0; ip < nParts; ip++)
{
    Info << "theta" << ip << " post Var transfer: [" << min(theta[ip]).value() << ", " << max(theta[ip]).value() << "]" << endl;
}
for(label ip = 0; ip < nParts; ip++)
{
    Info << "thetaVar" << ip << " post Var transfer: [" << min(thetaVar[ip]).value() << ", " << max(thetaVar[ip]).value() << "]" << endl;
}
Info << "TOTAL INTERNAL ENERGY POST TRANSFER: " << sum(theta.sum()).value() << endl;
