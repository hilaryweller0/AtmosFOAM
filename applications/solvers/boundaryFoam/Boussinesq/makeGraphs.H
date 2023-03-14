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
    turbulence->sigma()
);

// Lookup Cmu corresponding to the turbulence model selected
const scalar Cmu = turbulence->coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);
//const scalar s = turbulence->coeffDict().lookupOrDefault<scalar>
//    ("kappa", 0.41);

volScalarField mixingLength
(
    "mixingLength",
    turbulence->nut()/(Cmu/*::pow(Cmu,0.25)*/ * sqrt(turbulence->k()))
);

runTime.write();

setWriter::New(runTime.graphFormat())->write
(
    runTime.globalPath()
   /functionObjects::writeFile::outputPrefix
   /args.executable()
   /runTime.timeName(),

    args.executable(),

    coordSet(true, "y", y),

    "U", (flowDirection & U)().primitiveField(),

    "nu", turbulence->nu()().primitiveField(),
    "k", turbulence->nu()().primitiveField(),
    "epsilon", turbulence->nu()().primitiveField(),

    "Rff", (flowDirection & R & flowDirection)().primitiveField(),
    "Rww", (wallNormal & R & wallNormal)().primitiveField(),
    "Rfw", (flowDirection & R & wallNormal)().primitiveField(),

    "u", sqrt(mag(R.component(symmTensor::XX)))().primitiveField(),
    "v", sqrt(mag(R.component(symmTensor::YY)))().primitiveField(),
    "w", sqrt(mag(R.component(symmTensor::ZZ)))().primitiveField(),
    "uy", R.component(symmTensor::XY)().primitiveField(),

    "gammaDot", mag(fvc::grad(U))().primitiveField(),
    "lm", mixingLength.primitiveField(),

    "TKEproduction", (magSqr(fvc::grad(U))*turbulence->nut())().primitiveField(),
    
    "TKEturbulentTransport",
    fvc::laplacian
    (
        turbulence->nut()
            /turbulence->coeffDict().lookupOrDefault<scalar>("ksigma",1),
        turbulence->k()
    )().primitiveField(),
    
    "TKEviscousTransport",
    fvc::laplacian(turbulence->nu(), turbulence->k())().primitiveField(),
    
    "b", b.primitiveField(),
    "heatFlux", (Prandtl*turbulence->nut()*fvc::grad(b))().primitiveField()
);

