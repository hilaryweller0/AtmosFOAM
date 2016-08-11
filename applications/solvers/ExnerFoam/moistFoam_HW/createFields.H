Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info << "Reading in air and water phases" << endl;
wordList phaseNames(2);
phaseNames[0] = "air";
phaseNames[1] = "water";
atmosphere atmos(phaseNames, runTime.timeName(), mesh, thermoDict);

perfectGasPhase& air = atmos[0].gas();
fluidSpecie& water = atmos[1];

Info << "Reading in theta" << endl;
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
theta.oldTime();

Info << "Calculating, density, Exner, pressure and temperature from atmos and theta" << endl;

volScalarField rho("rho", atmos.sumDensity());
volScalarField Exner("Exner", atmos.ExnerFromTheta(theta));
volScalarField T("T", Exner*theta);
volScalarField p("p", air.pFromExner(Exner));

rho.write();
Exner.write();
T.write();
p.write();

Info << "Reading in velocity field" << endl;
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

// Fields needed in the Helmholtz equation
volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

// fixed flux boundary conditions
wordList VuBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "fixedValue";
    }
}

// Cell centre velocity
volVectorField u
(
    IOobject("u", runTime.timeName(), mesh),
    H.reconstruct(Uf & mesh.Sf()),
    uBCs
);

// Velocity flux
surfaceScalarField un("un", Uf & mesh.Sf());
volScalarField divu("divu", fvc::div(un));

