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

surfaceScalarField volFlux
(
    IOobject
    (
        "volFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Uf & mesh.Sf()
);

surfaceScalarField massFlux("massFlux", volFlux*linearInterpolate(T[0]));
surfaceScalarField* phip = &volFlux;
#define phi (*phip)

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    CourantNo(volFlux, runTime.deltaT())
);
Co.write();

surfaceScalarField beta("beta", max(scalar(0), 1-CoLimit/maxInterp.interpolate(Co)));
surfaceScalarField alpha("alpha", max(scalar(0.5), beta));
Info << "alpha goes from " << min(alpha).value() << " to "
     << max(alpha).value() << endl;
Info << "beta goes from " << min(beta).value() << " to "
     << max(beta).value() << endl;

dimensionedScalar Vtot = sum(mesh.V());

// Additional fluxes for each tracer for RK time stepping with FCT
PtrList<surfaceScalarField> fluxOld(T.size());
PtrList<surfaceScalarField> fluxT(T.size());
PtrList<PtrList<surfaceScalarField>> fluxTmp(T.size());

for(label iT = 0; iT < T.size(); iT++)
{
    fluxOld.set(iT, 0*phi*linearInterpolate(T[iT]));
    fluxT.set(iT, fluxOld[iT]);
    fluxTmp.set(iT, new PtrList<surfaceScalarField>(Bt.nSteps()));
    for (int iRK = 0; iRK < Bt.nSteps(); iRK++)
    {
        fluxTmp[iT].set(iRK, fluxOld[iT]);
    }
    if (withDensity) phip = &massFlux;
}
if (withDensity) phip = &volFlux;

