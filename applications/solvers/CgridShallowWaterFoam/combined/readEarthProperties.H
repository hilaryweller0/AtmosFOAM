Info<< "\nReading earthProperties" << endl;

IOdictionary earthProps
(
    IOobject
    (
        "earthProperties",
         runTime.constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
    )
);

const dimensionedScalar magg(earthProps.lookup("magg"));
const Switch rotating(earthProps.lookup("rotating"));
const dimensionedVector Omega =
    rotating
  ? earthProps.lookup("Omega")
  : dimensionedVector("Omega", dimless/dimTime, vector::zero);

Info<< "Reading field h0 (mountain heights) if present\n" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
         runTime.constant(),
         mesh,
         IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar(dimLength, 0)
);
Info << "Mountain height from " << min(h0.primitiveField()) << " to "
     << max(h0.primitiveField()) << endl;

const volVectorField gHat = mesh.C()/mag(mesh.C());
const surfaceVectorField gHatf = mesh.Cf()/mag(mesh.Cf());
     
Info<< "Creating Coriolis Force" << endl;
const surfaceVectorField F("F", ((2.0*Omega) & gHatf)*gHatf);

// Mean depth for linearised SWE
const dimensionedScalar H(earthProps.lookup("H"));
