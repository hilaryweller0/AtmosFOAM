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

/*volVectorField U
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
*/
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh //Uf & mesh.Sf()
);

/*volScalarField rho
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
surfaceScalarField phiv("phiv", phi/rhof);*/
//surfaceScalarField& phiv(phi);

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    /*CourantNo(phi/rhof, runTime.deltaT())*/
    CourantNo(phi, runTime.deltaT())
);
Co.write();

// Courant number on faces
surfaceScalarField Cof("Cof", maxInterp.interpolate(Co));

// Implicit/explicit dependent on Courant number
surfaceScalarField ImEx("ImEx", 0.5*(sign(Cof - CoLimit) + 1));

// Temporal off centering that varies in space (large for large Courant)
surfaceScalarField offCentre("offCentre", max(0.5, 1 - 1/Cof));

