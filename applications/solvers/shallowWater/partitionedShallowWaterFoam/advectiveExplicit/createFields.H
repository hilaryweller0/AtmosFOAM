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

//for (label ip=0; ip<nParts; ip++)
//{
//    u[ip] = fvc::reconstruct(volFlux[ip]);
//}

partitionedVolVectorField u_temp("u_temp", partNames, u[0], sigma);
partitionedVolVectorField hu("hu", partNames, sigmah[0]*u[0]);
partitionedVolVectorField hu_temp("hu_temp", partNames, hu[0]);

for (label ip = 1; ip < nParts; ip++)
{
    u_temp[ip] = u[ip];
    hu[ip] = sigmah[ip]*u[ip];
    hu_temp[ip] = hu[ip];
}

Info << "Initialising drag for each fluid." << endl;
partitionedSurfaceVectorField drag
(
    "drag", partNames, 0*Uf[0]/dt, sigmaf
);

partitionedVolScalarField sink
(
    "sink", partNames, 0*sigma[0]/dt, sigma
);
partitionedVolScalarField totalSink
(
    "totalSink", partNames, 0*sigma[0]/dt, sigma
);
if (useMassTransfer)
{
    Info << "Setting transfer terms." << endl;
    partitionedVolScalarField sink
    (
        "sink", partNames, mesh, runTime.timeName()
    );
    
    partitionedVolScalarField momentumSink
    (
        "momentumSink", partNames, mesh, runTime.timeName()
    );
}

partitionedVolScalarField momentumSource
(
    "momentumSource", partNames, 0*sigma[0]/dt, sigma
);
if (useBuoyancy)
{
    partitionedVolScalarField momentumSource
    (
        "momentumSource", partNames, mesh, runTime.timeName()
    );
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
    0.5*sigmah[0]*magSqr(u[0]) + 0.5*sigmah[1]*magSqr(u[1]) + 0.5*gravity*sqr(h)
);

// Store old time levels for some partitionedFields
sigmah.storeTime();
volFlux.storeTime();
Uf.storeTime();

sigmah.write();
energy.write();
