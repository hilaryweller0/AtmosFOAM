// So postProcess can run properly
#include "readEnvironmentalProps.H"
#define dt runTime.deltaT()

// Read in prognostic variables
Info<< "Reading field P\n" << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading in buoyancy, b\n" << endl;
volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading in cell centre velocity, u\n" << endl;
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

// fixed flux boundary conditions
wordList fluxBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = u.boundaryField()[patchi].type();
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

Info << "Calculating the volume flux\n"<<endl;
surfaceScalarField volFlux
(
    IOobject
    (
       "volFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(u) & mesh.Sf(),
    fluxBCs
);

// Velocity on the faces
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(u)
);
Uf.write();

// Surface scalar fields for hydrostatic boundary conditions
const surfaceScalarField gradPcoeff
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("1", dimless, scalar(1))
);

surfaceScalarField bf("bf", fvc::interpolate(b, "b")*mesh.Sf().component(2));

/*// Viscosity that may vary in space
const volScalarField nu
(
    IOobject("nu", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "nu",
        dimensionedScalar("nu", dimensionSet(0,2,-1,0,0), scalar(0))
    )
);
const dimensionedScalar maxNu = max(nu);*/

// Rate of change of b (eg for radiative cooling)
volScalarField Q
(
    IOobject("Q", runTime.timeName(), mesh, 
             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Q", dimensionSet(0,1,-3,0,0), scalar(0))
);

// Total initial heat in the domain
dimensionedScalar totalHeat0
(
    "totalHeat0", dimensionSet(0,1,-2,0,0), scalar(0)
);
totalHeat0 = fvc::domainIntegrate(b)/sum(mesh.V());

// Store old times and rates of change
Info << "Storing necessary old times and rates of change\n" << endl;
volFlux.oldTime();
b.oldTime();

volScalarField dbdt
(
    "dbdt",
   -fvc::div(volFlux, b, "b") + b*fvc::div(volFlux)
);
dbdt.oldTime();

surfaceScalarField dVolFluxDt
(
    IOobject
    (
        "dVolFluxDt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    nu*
    (
        (fvc::interpolate(fvc::laplacian(u)) & mesh.Sf())
      + fvc::snGrad(fvc::div(volFlux - (Uf & mesh.Sf())))*mesh.magSf()
    )
    -(!args.options().found("linear"))*
    (
        mesh.Sf() & (fvc::interpolate(fvc::div(volFlux, u, "u"), "convection"))
    )
  + bf
  - fvc::snGrad(P)*mesh.magSf()
  - (mesh.Sf() & (2*(Omega ^ Uf)))
);
dVolFluxDt.oldTime();

