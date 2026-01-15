Info<< "\nReading earthProperties" << endl;

IOdictionary earthProperties
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

const dimensionedScalar magg(earthProperties.lookup("magg"));
const Switch rotating(earthProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating
  ? earthProperties.lookup("Omega")
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