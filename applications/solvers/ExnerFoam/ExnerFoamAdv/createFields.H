// So postProcess can run properly
#include "readEnvironmentalProperties.H"
#include "readThermoProperties.H"
HodgeOps H(mesh);
#define dt runTime.deltaT()

// Read in prognostic variables
Info<< "Reading field Exner\n" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ExnerPrime
(
    IOobject("ExnerPrime", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("zero", dimless, scalar(0)),
    "zeroGradient"
);

Info << "Reading in theta\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
surfaceScalarField thetaf("thetaf", fvc::interpolate(theta));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);

Info << "Reading in the velocity\n" << endl;
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// fixed flux boundary conditions
wordList fluxBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        fluxBCs[patchi] = "fixedValue";
        uBCs[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        fluxBCs[patchi] = "fixedValue";
        uBCs[patchi] = "fixedValue";
    }
}

Info << "Calculating the mass flux and volume flux\n"<<endl;
surfaceScalarField volFlux
(
    IOobject
    (
       "volFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Uf & mesh.Sf(),
    fluxBCs
);
surfaceScalarField flux
(
    IOobject
    (
        "flux", 
        runTime.timeName(), 
        mesh,
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho)*volFlux
);

// Cell centre velocity
Info << "Reading in BCs and calculating the cell centre velocity, u\n" << endl;
volVectorField u
(
    IOobject
    (
        "u", 
        runTime.timeName(), 
        mesh,
        IOobject::MUST_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh
);
u = fvc::reconstruct(volFlux);

Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.oldTime();
flux.oldTime();
theta.oldTime();
rho.oldTime();

volScalarField dThetaDt
(
    "dThetaDt",
   -fvc::div(volFlux, theta) + theta*fvc::div(volFlux)
);
dThetaDt.oldTime();

surfaceScalarField dVolFluxDt("dVolFluxDt", volFlux/dt);
dVolFluxDt.oldTime();

