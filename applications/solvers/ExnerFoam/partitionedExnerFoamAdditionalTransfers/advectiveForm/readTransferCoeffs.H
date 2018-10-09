Info << "\nReading environmentalProperties" << endl;

IOdictionary transferProps
(
    IOobject
    (
        "transferProperties",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar cloudRadiusMax(transferProps.lookup("cloudRadiusMax"));
const dimensionedScalar cloudRadiusMin(transferProps.lookup("cloudRadiusMin"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));

// Diffusion coefficient for inter-partition mass transfer
const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

// Is mass transferred based on partition divergence
const Switch divTransfer(transferProps.lookup("divTransfer"));

// Is mass transfer decided based on vertical velocity
const Switch wTransfer(transferProps.lookup("wTransfer"));

// Is mass transfer decided based on grad(theta)
const Switch thetaTransfer(transferProps.lookup("thetaTransfer"));
const dimensionedScalar thetaTransferDiffusivity
(
    transferProps.lookup("thetaTransferDiffusivity")
);

// Minimum sigmaRho for each partition
const scalar minSigma = readScalar(transferProps.lookup("minSigma"));

//Numerical treatment of transfer terms
const scalar alphaM = readScalar
(
    transferProps.lookup("massTransferImplicit")
);
const scalar alphaV = readScalar
(
    transferProps.lookup("velocityTransferImplicit")
);
const scalar alphaT = readScalar
(
    transferProps.lookup("thetaTransferImplicit")
);
