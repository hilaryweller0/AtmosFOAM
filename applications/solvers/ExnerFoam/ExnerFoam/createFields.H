// Read in prognostic variables
Info<< "Reading field Uf\n" << endl;
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

// Cell centre velocity
//volVectorField U("U", H.reconstruct(Uf & mesh.Sf()));
volVectorField U("U", fvc::reconstruct(Uf & mesh.Sf()));

Info<< "Reading field Exner\n" << endl;
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

Info<< "Reading field theta\n" << endl;
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

// Old part of theta change (and theta on the face)
surfaceScalarField thetaf
(
    IOobject("thetaf", runTime.timeName(), mesh),
    fvc::interpolate(theta)
);
thetaf.oldTime();

// Reading in the sponge layer coefficients (if present)
surfaceScalarField muSponge
(
    surfaceScalarField
    (
        IOobject("muSponge", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    )/dt
);
Info << "Maximum sponge value = " << max(muSponge).value()*dt.value()
     << endl;

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

// fixed flux boundary conditions
wordList phiUBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || utype == "fixedValue"
    )
    {
        phiUBCs[patchi] = "fixedValue";
    }
}

Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
       "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhof*(Uf & mesh.Sf())
//    phiUBCs
);
phi.oldTime();

// Temporary flux with fixed boundary conditions (without hydrostatic terms)
surfaceScalarField phiU
(
    IOobject("phiU", runTime.timeName(), mesh),
    phi,
    phiUBCs
);

// Perpendicular component of the momentum
surfaceVectorField vperp("vperp", Uf - (Uf & mesh.Sf())*mesh.Sf()/sqr(mesh.magSf()));

// Rates of change for Crank-Nicholson
surfaceScalarField dPhidt
(
    IOobject("dPhidt", runTime.timeName(), mesh),
    -(mesh.Sf() & (fvc::interpolate(fvc::div(phi,U))+2*rhof*(Omega^(Uf-Ug))))
    +rhof*(gSf - mesh.magSf()*Cp*thetaf*fvc::snGrad(Exner)),
    phiUBCs
);

volScalarField divPhi("divPhi", fvc::div(phi));
volScalarField divPhitheta("divPhitheta", fvc::div(phi, theta));
divPhi.oldTime();
divPhitheta.oldTime();

// fix no flow boundaries
forAll(phiU.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = U.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
    )
    {
        phiU.boundaryFieldRef()[patchi] = 0;
        dPhidt.boundaryFieldRef()[patchi] = 0;
    }
}

