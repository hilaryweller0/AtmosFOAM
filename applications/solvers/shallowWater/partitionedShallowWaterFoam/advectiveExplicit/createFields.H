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

volScalarField gravity
(
    IOobject
    (
        "gravity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Setting up fluid types." << endl;
const label nParts = fluidsDict.lookupOrDefault<label>
(
    "numberOfFluids", 2
);
wordList partNames(nParts);
std::ostringstream fluidDefaultName;
std::ostringstream fluidSettingName;

Info << "Generating fields for " << nParts << " fluids" << endl;
for (label ip=0; ip<nParts; ip++)
{
    fluidSettingName.str(""); 
    fluidSettingName.clear();
    fluidSettingName << "fluidName" << ip;
    
    fluidDefaultName.str(""); 
    fluidDefaultName.clear();
    fluidDefaultName << "fluid" << ip;
    
    partNames[ip] = fluidsDict.lookupOrDefault<string>
    (
        fluidSettingName.str(), fluidDefaultName.str()
    );
    Info << fluidSettingName.str() << " " << fluidDefaultName.str() << endl;
    Info << "Fluid " << ip << " name: " << partNames[ip] << endl;
}

Info << "Reading sigma (volume fraction) for each fluid." << endl;
partitionedVolScalarField sigma
(
    "sigma", partNames, mesh, runTime.timeName()
);

partitionedSurfaceScalarField sigmaf
(
    "sigmaf", partNames, fvc::interpolate(sigma[0], "sigma")
);
for (label ip=0; ip<nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
}
sigmaf.updateSum();

Info << "Setting effective fluid heights." << endl;
partitionedVolScalarField sigmah
(
    "sigmah", partNames, sigma[0] * h
);
sigmah.updateSum();
partitionedSurfaceScalarField sigmahf
(
    "sigmahf", partNames, fvc::interpolate(sigmah[0])
);
partitionedVolScalarField sigmah_temp
(
    "sigmah_temp", partNames, sigmah[0]
);
for (label ip=0; ip<nParts; ip++)
{
    sigmah[ip] = sigma[ip]*h;
    sigmahf[ip] = fvc::interpolate(sigmah[ip],"sigmah");
    sigmah_temp[ip] = sigmah[ip];
}

Info << "Reading velocity of each fluid." << endl;
partitionedSurfaceVectorField Uf
(
    "Uf", partNames, mesh, runTime.timeName(), sigmaf
);

Info << "Diagnosing fluxes for each fluid." << endl;
partitionedSurfaceScalarField volFlux
(
    //"volFlux", partNames, Uf[0] & mesh.Sf(), sigmaf
    "volFlux", partNames, mesh, runTime.timeName(), sigmaf
);
partitionedSurfaceScalarField flux
(
    "flux", partNames, fvc::interpolate(sigmah[0])*volFlux[0]
);
for (label ip=0; ip<nParts; ip++)
{
    volFlux[ip] = Uf[ip] & mesh.Sf();
    flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
}
volFlux.updateSum();

Info << "Diagnosing cell-centred velocity of each fluid." << endl;
partitionedVolVectorField u
(
    //"u", partNames, fvc::reconstruct(volFlux[0]), sigma
    "u", partNames, mesh, runTime.timeName(), sigma
);
for (label ip=0; ip<nParts; ip++)
{
    u[ip] = fvc::reconstruct(volFlux[0]);
}

Info << "Initialising drag for each fluid." << endl;
partitionedSurfaceVectorField drag
(
    "drag", partNames, 0*Uf[0]/dt, sigmaf
);

Info << "Setting transfer terms." << endl;
partitionedVolScalarField sink
(
    "sink", partNames, mesh, runTime.timeName()
);
partitionedVolScalarField totalSink
(
    "totalSink", partNames, 0*sigmah[0], sigma
);
partitionedVolScalarField momentumSink
(
    "momentumSink", partNames, mesh, runTime.timeName()
);
partitionedVolScalarField momentumSource
(
    "momentumSource", partNames, mesh, runTime.timeName()
);
