Info<< "Reading tracers\n" << endl;

PtrList<volScalarField> T(Tnames.size());
for(label iT = 0; iT < T.size(); iT++)
{
    Info << "Reading " << Tnames[iT] << endl;
    T.set
    (
        iT, 
        volScalarField
        (
            IOobject(Tnames[iT], runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        )
    );

    T[iT].writeOpt() = IOobject::AUTO_WRITE;
    T[iT].checkIn();
    T[iT].rename(Tnames[iT]);
    T[iT].correctBoundaryConditions();
    T[iT].oldTime();
}

volScalarField Tsum = withDensity ?
    volScalarField
    (
        IOobject(T[0].name()+T[1].name(), runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        T[0]*T[1]
    ) :
    volScalarField
    (
        IOobject("Tsum", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        T[0]
    );
if (!withDensity) for(label iT = 1; iT < T.size(); iT++)
{
    Tsum += T[iT];
}
Tsum.write();

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
phi.oldTime();

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(phi, runTime.deltaT())
);
Co.write();

surfaceScalarField beta("beta", max(scalar(0), 1-CoLimit/maxInterp.interpolate(Co)));
surfaceScalarField alpha("alpha", max(scalar(0.5), beta));

dimensionedScalar Vtot = sum(mesh.V());

PtrList<volScalarField> divSum(Tnames.size());
for(label iT = 0; iT < divSum.size(); iT++)
{
    divSum.set(iT, fvc::div(phi*linearInterpolate(T[iT])));
}
