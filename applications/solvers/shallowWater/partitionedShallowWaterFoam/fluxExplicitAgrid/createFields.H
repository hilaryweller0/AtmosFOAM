// Read in prognostic variables
Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
h.oldTime();

// The partition names
const label nParts = 2;
wordList partName(nParts);
partName[0] = "stable";
partName[1] = "buoyant";
wordList partNames(nParts);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<volScalarField> sigmah(nParts);
PtrList<volScalarField> sigmah_temp(nParts);
PtrList<volScalarField> hOld(nParts);
PtrList<surfaceScalarField> sigmahf(nParts);
PtrList<surfaceVectorField> Uf(nParts);
PtrList<volVectorField> u(nParts);
PtrList<volVectorField> hu(nParts);
PtrList<volVectorField> hu_temp(nParts);
PtrList<surfaceScalarField> flux(nParts);
PtrList<surfaceScalarField> volFlux(nParts);
PtrList<volScalarField> sink(nParts);
PtrList<volScalarField> totalSink(nParts);
PtrList<volScalarField> momentumSource(nParts);

// Read in or set variables in partitions
for (label ip = 0; ip < nParts; ip++)
{
    sigma.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigma",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    
    sigmah.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigmah",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            h*sigma[ip]
        )
    );
    sigmah[ip].oldTime();

    sigmah_temp.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigmah_temp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            h*sigma[ip]
        )
    );

    sigmahf.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"sigmahf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            linearInterpolate(sigmah[ip])
        )
    );
    
    Uf.set
    (
        ip,
        new surfaceVectorField
        (
            IOobject
            (
                partNames[ip]+"Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    flux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"flux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            (Uf[ip] & mesh.Sf())*sigmahf[ip]
        )
    );
    flux[ip].oldTime();
    
    volFlux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"volFlux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Uf[ip] & mesh.Sf()
        )
    );

    u.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"u",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::reconstruct(flux[ip]/sigmahf[ip]),
            Uf[ip].boundaryField().types()
        )
    );
    
    hu.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"hu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigmah[ip]*u[ip]
        )
    );
    
    hu_temp.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"hu_temp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigmah[ip]*u[ip]
        )
    );
    
    sink.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sink",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
    
    totalSink.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"totalSink",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            0*sigma[0]/dt
        )
    );
    
    momentumSource.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"momentumSource",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
}

volScalarField delta_ke_centre
(
    IOobject
    (
        "delta_ke_centre",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.*sigmah[0]*magSqr(u[0])
);
