Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
T.correctBoundaryConditions();

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector::zero),
    "zeroGradient"
);

surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Uf & mesh.Sf()
);

Info << "Reading in density if present" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("rho", dimless, scalar(1))
);
surfaceScalarField rhof
(
    IOobject("rhof", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    linearInterpolate(rho)
);
// Volume flux (for Courant number limit)
surfaceScalarField phiv("phiv", phi/rhof);

// Check that phi dimensions are density * volFlux
if (phiv.dimensions() != dimensionSet(0,3,-1,0,0))
{
    FatalErrorIn("implicitAdvectionFoam")
        << "Dimension inconsistency between flux, phi = " << phi.dimensions()
        << "\nand density, rho = " << rho.dimensions() << exit(FatalError);
}

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(phiv, runTime.deltaT())
);
Co.write();
dimensionedScalar Vtot = sum(mesh.V());

