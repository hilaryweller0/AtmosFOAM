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

volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("rho", dimless, scalar(1))
);
surfaceScalarField rhof("rhof", linearInterpolate(rho));

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(phi, runTime.deltaT())
);
Co.write();

// Split the flux into parts for implicit and explicit solution
surfaceScalarField phiEx = 0*phi;

const surfaceScalarField phiLimit = 0.5*CoLimit*mesh.magSf()/mesh.deltaCoeffs()/dt;
if (CoLimit > SMALL) // Default is implicit with CoLimit = 0
{
   forAll(phi, faceI)
    {
        if (mag(phi[faceI]) < phiLimit[faceI])
        {
            phiEx[faceI] = phi[faceI];
            phi[faceI] = 0;
        }
    }
}

volScalarField CoEx
(
    IOobject("CoEx", runTime.timeName(), mesh),
    CourantNo(phiEx, dt)
);
CoEx.write();
volScalarField CoIm
(
    IOobject("CoIm", runTime.timeName(), mesh),
    CourantNo(phi, dt)
);
CoIm.write();

