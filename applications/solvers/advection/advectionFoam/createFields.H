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
// Ensure that T.oldTime() and T.oldTime.oldTime() are stored
T.oldTime();
T.oldTime().oldTime();

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
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

volScalarField k1
(
    IOobject
    (
        "k1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T)
);

volScalarField k2
(
    IOobject
    (
        "k2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T)
);

volScalarField k3
(
    IOobject
    (
        "k3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T)
);

volScalarField k4
(
    IOobject
    (
        "k4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phi,T)
);

// Calculate divU for checking divergence free
volScalarField divU("divU", fvc::div(phi));
divU.write();

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh),
    CourantNo(phi, dt)
);
Co.write();
if (args.options().found("implicitWhereNeeded"))
{
    // Split the flux into the large part and the smaller part
    surfaceScalarField phiSmall 
        = min(phi, 0.5*mesh.magSf()/mesh.deltaCoeffs()/dt);
    phiSmall = max(phiSmall, -0.5*mesh.magSf()/mesh.deltaCoeffs()/dt);
    surfaceScalarField phiBig = phi - phiSmall;

    volScalarField CoExp
    (
        IOobject("CoExp", runTime.timeName(), mesh),
        CourantNo(phiSmall, dt)
    );
    CoExp.write();
    volScalarField CoImp
    (
        IOobject("CoImp", runTime.timeName(), mesh),
        CourantNo(phiBig, dt)
    );
    CoImp.write();
}
