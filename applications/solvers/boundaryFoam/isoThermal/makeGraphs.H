volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->R()
);

// Lookup Cmu corresponding to the turbulence model selected
const scalar Cmu = turbulence->coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);
const scalar kappa = turbulence->coeffDict().lookupOrDefault<scalar>
    ("kappa", 0.41);

volScalarField mixingLength
(
    "mixingLength",
    turbulence->nut()/(::pow(Cmu,0.25) * sqrt(turbulence->k()))
);

wallDist wd(mesh);
volScalarField mixingLength2("mixingLength2", kappa*wd.y());

runTime.write();

const word& gFormat = runTime.graphFormat();

makeGraph(y, flowDirection & U, "Uf", gFormat);

makeGraph(y, turbulence->nut(), gFormat);
makeGraph(y, turbulence->k(), gFormat);
makeGraph(y, turbulence->epsilon(), gFormat);
//makeGraph(y, turbulence->omega(), gFormat);

makeGraph(y, flowDirection & R & flowDirection, "Rff", gFormat);
makeGraph(y, wallNormal & R & wallNormal, "Rww", gFormat);
makeGraph(y, flowDirection & R & wallNormal, "Rfw", gFormat);

makeGraph(y, sqrt(mag(R.component(symmTensor::XX))), "u", gFormat);
makeGraph(y, sqrt(mag(R.component(symmTensor::YY))), "v", gFormat);
makeGraph(y, sqrt(mag(R.component(symmTensor::ZZ))), "w", gFormat);
makeGraph(y, R.component(symmTensor::XY), "uv", gFormat);

makeGraph(y, mag(fvc::grad(U)), "gammaDot", gFormat);
makeGraph(y, mixingLength, "lm", gFormat);
makeGraph(y, mixingLength2, "lm2", gFormat);
makeGraph(y, sqr(mixingLength2)*mag(fvc::grad(U)), "nut2", gFormat);

makeGraph(y, magSqr(fvc::grad(U))*turbulence->nut(), "TKEproduction", gFormat);
makeGraph
(
    y,
    fvc::laplacian
    (
        turbulence->nut()/turbulence->coeffDict().lookupOrDefault<scalar>("ksigma",1),
        turbulence->k()
    ),
    "TKEturbulentTransport",
    gFormat
);

makeGraph
(
    y,
    fvc::laplacian(turbulence->nu(), turbulence->k()),
    "TKEviscousTransport",
    gFormat
);

