Info << "Reading in Exner pressure" << endl;
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
Exner.oldTime();

Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info << "Reading in air and water phases in stable and buoyant partitions"
     << endl;
wordList phaseNames(2);
phaseNames[0] = "air";
phaseNames[1] = "water";
wordList partNames(2);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

partitionedAtmosphere atmosParts(partNames, phaseNames, mesh, thermoDict,Exner);

// Select just the stable part
partition& atmos = atmosParts[0];

perfectGasPhase& air = atmos[0].gas();
//fluidSpecie& water = atmos[1];


Info << "Calculating the total density, pressure, temperature and velocity"
     << endl;

//volScalarField rho("rho", atmosParts.rho());
//volScalarField theta("theta", atmosParts.rhoTheta()/rho);
//volScalarField T("T", theta*Exner);
volScalarField p("p", air.pFromExner(Exner));
//surfaceVectorField Uf("Uf", atmosParts.rhoUf()/linearInterpolate(rho));

//// Time-stepping variables
//rho.oldTime();
//volScalarField dRhodt("dRhodt", 0*rho/dt);
//dRhodt.oldTime();

//volScalarField dRhoThetadt("dRhoThetadt", 0*rho*theta/dt);
//dRhoThetadt.oldTime();
//volScalarField rhoThetaSource("rhoThetaSource", dRhoThetadt);
//surfaceScalarField rhof("rhof", fvc::interpolate(rho));


// Fields needed in the Helmholtz equation
//volScalarField Psi("Psi", rho/Exner);
//Psi.oldTime();
//surfaceScalarField gradPcoeff
//(
//    "gradPcoeff",
//    fvc::interpolate
//    (
//        air.Cp()*theta*atmos.rhoR()
//       /(rho*air.R()*atmos.volGas())
//    )
//);

//// fixed flux boundary conditions
//wordList VuBCs(mesh.boundaryMesh().size(), "calculated");
//wordList uBCs(mesh.boundaryMesh().size(), "zeroGradient");
//forAll(mesh.boundaryMesh(), patchi)
//{
//    const word btype = mesh.boundaryMesh()[patchi].type();
//    const word utype = Uf.boundaryField()[patchi].type();
//    if (btype == "wall" || btype == "symmetryPlane" || btype == "empty")
//    {
//        VuBCs[patchi] = "fixedValue";
//        uBCs[patchi] = "slip";
//    }
//    else if (utype == "fixedValue")
//    {
//        VuBCs[patchi] = "fixedValue";
//        uBCs[patchi] = "fixedValue";
//    }
//}

//// Cell centre velocity
//volVectorField u
//(
//    IOobject("u", runTime.timeName(), mesh),
//    H.reconstruct(Uf & mesh.Sf()),
//    uBCs
//);

//Info<< "Calculating face flux field U\n" << endl;
//surfaceScalarField U
//(
//    IOobject
//    (
//       "U",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    rhof*(Uf & mesh.Sf())
////    phiUBCs
//);
//U.oldTime();

//Info << "Calculating V mass circulation field\n" << endl;
//surfaceScalarField V
//(
//    IOobject("V", runTime.timeName(), mesh),
//    rhof*(Uf & H.delta())
////    phiUBCs
//);
//V.oldTime();

//// Velocity flux
//volScalarField divu("divu", fvc::div(Uf & mesh.Sf()));

//// Temporary flux with fixed boundary conditions (without hydrostatic terms)
//surfaceScalarField Vu
//(
//    IOobject("Vu", runTime.timeName(), mesh),
//    V,
//    VuBCs
//);

//// Rates of change for Crank-Nicholson
//surfaceScalarField dVdt
//(
//    IOobject("dVdt", runTime.timeName(), mesh),
//    -(H.delta() & (fvc::interpolate(fvc::div(U,u))+2*rhof*(Omega^(Uf-Ug))))
//    +rhof*(gd - H.magd()*gradPcoeff*fvc::snGrad(Exner)),
//     VuBCs
//);
//dVdt.oldTime();

//// Reading in the sponge layer coefficients (if present)
//surfaceScalarField muSponge
//(
//    surfaceScalarField
//    (
//        IOobject("muSponge", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
//        mesh,
//        dimensionedScalar("muSponge", dimless, scalar(0))
//    )/dt
//);
//Info << "Maximum sponge value = " << max(muSponge).value()*dt.value()
//     << endl;

