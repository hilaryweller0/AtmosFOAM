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

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh),
    CourantNo(phi, dt)
);
Co.write();

surfaceScalarField phiSmall = phi;
surfaceScalarField phiBig = 0*phi;

// Split the flux into the large part and the smaller part
surfaceScalarField phiLimit = 0.5*CoLimit*mesh.magSf()/mesh.deltaCoeffs()/dt;
forAll(phi, faceI)
{
    if (mag(phi[faceI]) > phiLimit[faceI])
    {
        phiSmall[faceI] = 0;
        phiBig[faceI] = phi[faceI];
    }
}

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

volScalarField divPhiOld = fvc::div(phiSmall, T, "explicit")
                         + fvc::div(phiBig, T, "implicit");

