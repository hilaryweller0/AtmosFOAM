// Read in prognostic variables
Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
h.oldTime();

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<surfaceVectorField> Uf(nParts);
PtrList<volVectorField> u(nParts);
PtrList<surfaceScalarField> flux(nParts);
PtrList<surfaceScalarField> dFluxdt(nParts);
PtrList<surfaceScalarField> sigmahf(nParts);

// Read in or set variables in partitions
for (label ip = 0; ip < nParts; ip++)
{
    sigma.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigma",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    sigma[ip].oldTime();

    Uf.set
    (
        ip,
        new surfaceVectorField
        (
            IOobject
            (
                partNames[ip]+"Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    u.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"u",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            H.reconstruct(Uf[ip] & mesh.Sf()),
            Uf[ip].boundaryField().types()
        )
    );
    
    sigmahf.set
    (
        ip,
        new surfaceScalarField 
        (
            IOobject
            (
                partNames[ip]+"sigmahf",
                runTime.timeName(),
                mesh
            ),
            linearInterpolate(sigma[ip]*h)
        )
    );
    
    flux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"flux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sigmahf[ip]*(Uf[ip] & mesh.Sf())
        )
    );
    sigmahf[ip] = fvc::interpolate(sigma[ip]*h);
    flux[ip] = sigmahf[ip]*(Uf[ip] & mesh.Sf());
    flux[ip].oldTime();
    
    dFluxdt.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"dFluxdt",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dFluxdt", dimensionSet(0,4,-2,0,0), scalar(0)),
            "fixedValue"
        )
    );
    dFluxdt[ip].oldTime();
}

surfaceScalarField fluxSum
(
    IOobject("fluxSum", runTime.timeName(), mesh),
    flux[0]
);
for(label ip = 1; ip < nParts; ip++){ fluxSum += flux[ip]; }

// Summed variables
volScalarField sumSigma
(
    IOobject
    (
        "sumSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimless, scalar(1))
);

// Find sigmahf which sums to hf
surfaceScalarField hf("hf", fvc::interpolate(h));
surfaceScalarField sumSigmahf("sumSigmahf", hf*0);
for(label ip = 0; ip < nParts; ip++)
{
    sigmahf[ip] = fvc::interpolate(sigma[ip]*h);
    sumSigmahf += sigmahf[ip];
}
for(label ip = 0; ip < nParts; ip++)
{
    sigmahf[ip] *= hf/sumSigmahf;
}

