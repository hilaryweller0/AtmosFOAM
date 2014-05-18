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

//// Diffusion coefficient for theta
//const surfaceScalarField diffuseTheta
//(
//    IOobject("diffuseTheta", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
//    mesh,
//    dimensionedScalar("diffuseTheta", dimensionSet(1,-1,-1,0,0), scalar(0))
//);

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
wordList VuBCs(mesh.boundaryMesh().size(), "calculated");
wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = Uf.boundaryField()[patchi].type();
    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "slip";
    }
    else if (utype == "fixedValue")
    {
        VuBCs[patchi] = "fixedValue";
        uBCs[patchi] = "fixedValue";
    }
}

// Cell centre velocity
volVectorField u
(
    IOobject("u", runTime.timeName(), mesh),
    H.reconstruct(Uf & mesh.Sf()),
    uBCs
);

Info<< "Calculating face flux field U\n" << endl;
surfaceScalarField U
(
    IOobject
    (
       "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhof*(Uf & mesh.Sf())
//    phiUBCs
);

Info << "Calculating V mass circulation field\n" << endl;
surfaceScalarField V
(
    IOobject("V", runTime.timeName(), mesh),
    rhof*(Uf & H.delta())
//    phiUBCs
);
V.oldTime();

// Temporary flux with fixed boundary conditions (without hydrostatic terms)
surfaceScalarField Vu
(
    IOobject("Vu", runTime.timeName(), mesh),
    V,
    VuBCs
);

// Perpendicular component of the momentum
surfaceVectorField vperp("vperp", rhof*(Uf - (Uf & H.delta())*H.delta()/sqr(H.magd())));

// Rates of change for Crank-Nicholson
surfaceScalarField dVdt
(
    IOobject("dVdt", runTime.timeName(), mesh),
    -(H.delta() & (fvc::interpolate(fvc::div(U,u))+2*rhof*(Omega^(Uf-Ug))))
    +rhof*gd - H.magd()*Cp*rhof*thetaf*fvc::snGrad(Exner),
     VuBCs
);
dVdt.oldTime();

volScalarField divU("divU", fvc::div(U));
volScalarField divUtheta("divUtheta", fvc::div(U, theta));
divU.oldTime();
divUtheta.oldTime();

// intermediate values for sub time-stepping of advection
surfaceScalarField divUu1
(
    IOobject("divUu1", runTime.timeName(), mesh),
    fvc::interpolate(fvc::div(U,u)) & H.delta(),
    VuBCs
);

// fix no flow boundaries
forAll(Vu.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word utype = U.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
    )
    {
        Vu.boundaryField()[patchi] == 0;
        dVdt.boundaryField()[patchi] == 0;
        divUu1.boundaryField()[patchi] == 0;
    }
}

surfaceScalarField divUu2("divUu2", divUu1);
surfaceScalarField divUu3("divUu3", divUu1);

//volScalarField divUt1("divUt1", fvc::div(U,theta));
//volScalarField divUt2("divUt2", divUt1);
//volScalarField divUt3("divUt3", divUt1);

