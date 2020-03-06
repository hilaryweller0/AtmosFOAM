Info<< "Reading field T\n" << endl;

Info << "Reading prognostic moisture variables. " << endl;
volScalarField rl
(
    IOobject
    (
        "rl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
rl.oldTime();

volScalarField rl_temp
(
    IOobject
    (
        "rl_temp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rl
);
rl_temp.oldTime();

volScalarField rv
(
    IOobject
    (
        "rv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
rv.oldTime();

volScalarField rv_temp
(
    IOobject
    (
        "rv_temp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rv
);
rv_temp.oldTime();

volScalarField rt
(
    IOobject
    (
        "rt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
rt.oldTime();

Info << "Reading diagnostic moisture variables. " << endl;
volScalarField rl_diag
(
    IOobject
    (
        "rl_diag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rv_diag
(
    IOobject
    (
        "rv_diag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rt_diag
(
    IOobject
    (
        "rt_diag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading analytic moisture variables. " << endl;
volScalarField rl_analytic
(
    IOobject
    (
        "rl_analytic",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rv_analytic
(
    IOobject
    (
        "rv_analytic",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rt_analytic
(
    IOobject
    (
        "rt_analytic",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading phase change variables. " << endl;
volScalarField Sl
(
    IOobject
    (
        "Sl",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rl*0/dt
);

volScalarField Sv
(
    IOobject
    (
        "Sv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rv*0/dt
);

Info << "Reading pressure field. " << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info << "Reading temperature field. " << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info << "Creating saturation vapour pressure field. " << endl;
const volScalarField es
(
    "es",
    es0*Foam::exp( 17.67*(T/T0-273.15)/(T/T0-29.65) )
);

Info << "Creating saturation vapour mixing ratio field. " << endl;
const volScalarField rvs
(
    "rvs",
    0.018/0.029 * es/(P-es)
);

Info << "Reading velocity fields. " << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

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
    fvc::reconstruct(phi)
);
U.write();

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
Uf.write();

