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

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(phi/rhof, runTime.deltaT())
);
Co.write();

// Split the flux into parts for implicit and explicit solution
const surfaceScalarField phiLimit = CoLimit*mesh.magSf()/mesh.deltaCoeffs()/dt;
surfaceScalarField exp = 0.5*(sign(phiLimit - mag(phiv)) + 1);
Info << "phiLimit goes from " << min(phiLimit.internalField()).value()
     << " to " << max(phiLimit.internalField()).value() << endl;
Info << "phiv goes from " << min(phiv.internalField()).value()
     << " to " << max(phiv.internalField()).value() << endl;
Info << "exp goes from " << min(exp.internalField()).value()
     << " to " << max(exp.internalField()).value() << endl;

