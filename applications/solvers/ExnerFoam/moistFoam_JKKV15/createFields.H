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

Info<< "Reading field thetaRho\n" << endl;
volScalarField thetaRho
(
    IOobject
    (
        "thetaRho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
thetaRho.oldTime();

// Water vapour and liquid water mass mixing ratios
volScalarField qv
(
    IOobject
    (
        "qv", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField ql
(
    IOobject
    (
        "ql", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);
ql.oldTime();
qv.oldTime();

// Related mixing ratios
volScalarField qd("qd", 1 - qv - ql);
volScalarField rv("rv", qv/qd);
volScalarField rl("rl", ql/qd);

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

// Thermodynamics variable kappa as a function of mixing ratios
volScalarField cpml("cpml", Cp + Cpv*rv + rl*Cpl);
volScalarField Rm("Rm", R + Rv*rv);
volScalarField kappam("kappam", Rm/cpml);

Info<< "Creating field rho from eqn of state\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    moistGasExnerRho(Exner, thetaRho, kappam, R, pRef)
);
rho.oldTime();
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

// Fields needed in the Helmholtz equation
volScalarField Psi("Psi", rho/Exner);
Psi.oldTime();

surfaceScalarField gradPcoeff = gradPCoeff(cpml, thetaRho, rv, epsilon);
surfaceScalarField gradPcoeff2("gradPcoeff2", gradPcoeff*rhof*dt);

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

// Velocity flux
//surfaceScalarField Uv("Uv", U/rhof);

// Perpendicular component of the momentum
//surfaceVectorField vperp("vperp", Uf - (Uf & H.delta())*H.delta()/sqr(H.magd()));

// Rates of change for Crank-Nicholson
surfaceScalarField dVdt
(
    IOobject("dVdt", runTime.timeName(), mesh),
    -(H.delta() & (fvc::interpolate(fvc::div(U,u))+2*rhof*(Omega^(Uf-Ug))))
    +rhof*
    (
        gd
      - H.magd()*gradPcoeff*
        (
            fvc::snGrad(Exner)
          - fvc::interpolate(Exner*Foam::log(Exner)/kappam)*fvc::snGrad(kappam)
        )
    ),
     VuBCs
);
dVdt.oldTime();

volScalarField drhodt("drhodt", -fvc::div(U));
drhodt.oldTime();
volScalarField dRhoThetadt("dRhoThetadt", -fvc::div(U, thetaRho));
dRhoThetadt.oldTime();
volScalarField condenseRate("condenseRate", qv/dt*0);
volScalarField drhovdt("drhovdt", -rho*condenseRate);
drhovdt.oldTime();
volScalarField drholdt("drhovdt", rho*condenseRate);
drholdt.oldTime();

// Primitive thermodynamic variables
volScalarField p("p", pFromExner(Exner, kappam, pRef));
volScalarField T("T", TFromThetaRho(thetaRho, Exner, rv, rl, epsilon));

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
        Vu.boundaryFieldRef()[patchi] == 0;
        dVdt.boundaryFieldRef()[patchi] == 0;
    }
}

// thetae for diagnostics
volScalarField Lv("Lv", latentHeat(T, Lv0, T0, Cpl, Cpv));
volScalarField thetae
(
    IOobject
    (
        "thetae", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    thetaeFromPrimitive(T, p, Lv, rv, rl, pRef, epsilon, Cp, Cpl, R)
);
