// Read in prognostic variablesInfo<< "\nReading hRef" << endl;
uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, 0)
);
Info<< "Calculating field g.h\n" << endl;
dimensionedScalar ghRef(- mag(g)*hRef);
volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

Info<< "Reading Exner and lnExner\n" << endl;
volScalarField lnExner
(
    IOobject
    (
        "lnExner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
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

Info<< "Reading theta\n" << endl;
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

Info<< "Reading T\n" << endl;
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

Info << "Reading/creating Exnerg\n" << endl;
volScalarField Exnerg
(
    IOobject
    (
        "Exnerg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    Exner - gh/(Cp*theta)
);

