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
    "sigmah", partNames, mesh, runTime.timeName()
);
sigmah.updateSum();

Info << "Diagnosing cell-centred velocity of each fluid." << endl;
partitionedVolVectorField u
(
    "u", partNames, mesh, runTime.timeName(), sigma
);

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
partitionedSurfaceScalarField volFlux_temp
(
    "volFlux_temp", partNames, volFlux[0]
);
partitionedSurfaceScalarField flux
(
    "flux", partNames, fvc::interpolate(sigmah[0])*volFlux[0]
);
partitionedSurfaceScalarField flux_temp
(
    "flux_temp", partNames, fvc::interpolate(sigmah[0])*volFlux[0]
);
for (label ip=0; ip<nParts; ip++)
{
    Uf[ip] = fvc::interpolate(u[ip]);
    volFlux[ip] = Uf[ip] & mesh.Sf();
    volFlux_temp[ip] = volFlux[ip];
    flux[ip] = fvc::interpolate(sigmah[ip])*volFlux[ip];
    flux_temp[ip] = flux[ip];
}
volFlux.updateSum();


partitionedVolScalarField SIGMAH
(
    "SIGMAH", partNames, mesh, runTime.timeName()
);
partitionedVolVectorField U
(
    "U", partNames, mesh, runTime.timeName(), sigma
);

Info << "Reading velocity of each fluid." << endl;
partitionedSurfaceVectorField UF
(
    "UF", partNames, mesh, runTime.timeName(), sigmaf
);

Info << "Diagnosing fluxes for each fluid." << endl;
partitionedSurfaceScalarField VOLFLUX
(
    "VOLFLUX", partNames, mesh, runTime.timeName(), sigmaf
);

for (label ip=0; ip<nParts; ip++)
{
    UF[ip] = fvc::interpolate(U[ip]);
    VOLFLUX[ip] = UF[ip] & mesh.Sf();
}

surfaceScalarField delta_ke_face
(
    IOobject
    (
        "delta_ke_face",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    0.*fvc::interpolate(sigmah[0])*volFlux[0]*volFlux[0]
);

volScalarField delta_ke_centre
(
    IOobject
    (
        "delta_ke_centre",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    0.*sigmah[0]*magSqr(u[0])
);

volScalarField energy
(
    IOobject
    (
        "energy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.5*gravity*sqr(sigmah[0]+sigmah[1])
);

for (label ip=0; ip<nParts; ip++)
{
    energy += sigmah[ip]*(U[ip] & u[ip]);
    energy += 0.5*SIGMAH[ip]*(u[ip] & u[ip]);
}
energy.write();

// Store old time levels for some partitionedFields
sigmah.storeTime();
volFlux.storeTime();
Uf.storeTime();
