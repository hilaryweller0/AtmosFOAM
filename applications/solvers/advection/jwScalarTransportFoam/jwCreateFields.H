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


Info<< "Reading field phi\n" << endl;

surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

// Read U if present, otherwise create and write (not used)
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
    fvc::reconstruct(phi) //Why do this?
);
U.write();

// Read Uf if present, otherwise create and write (not used)
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
    linearInterpolate(fvc::reconstruct(phi)) //why reconstruct from phi and what does this do?
);
Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf()); //derivation?
Uf.write();    //where does write to

// Store div(phi,T) for time-stepping
volScalarField divPhiT = fvc::div(phi, T); // is this ok with (phi,T,"lowOrder")
volScalarField divPhiTlo = fvc::div(phi, T, "lowOrder");
volScalarField divPhiThi = fvc::div(phi, T, "highOrder"); //the aim is to not use this
volScalarField divPhiTcorr = fvc::div(phi, T, "correction");
T.oldTime();
T.oldTime().oldTime();
divPhiT.oldTime();
divPhiT.oldTime().oldTime();
divPhiThi.oldTime();
divPhiThi.oldTime().oldTime();
divPhiTlo.oldTime();
divPhiTlo.oldTime().oldTime();
divPhiTcorr.oldTime();
divPhiTcorr.oldTime().oldTime();




