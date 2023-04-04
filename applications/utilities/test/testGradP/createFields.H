Info<< "Reading P\n" << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
P.correctBoundaryConditions();

Info<< "Reading Pf\n" << endl;
surfaceScalarField Pf
(
    IOobject
    (
        "Pf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(P)
);

IOdictionary tracerDict
(
    IOobject("tracerFieldDict", runTime.constant(), mesh, IOobject::READ_IF_PRESENT)
);
const dimensionedVector gradPexact
(
    "gradPexact",
    P.dimensions()/ dimLength,
    vector
    (
        readScalar(tracerDict.lookup("a_x")),
        readScalar(tracerDict.lookup("a_y")),
        readScalar(tracerDict.lookup("a_z"))
    )
);
/*IOdictionary tracerDict
(
    IOobject("expHeightTracerDict", runTime.constant(), mesh, IOobject::READ_IF_PRESENT)
);
volVectorField gradPexact
(
    IOobject("gradPexact", runTime.timeName(), mesh),
    mesh,
    dimensionedVector("gradPexact", P.dimensions()/ dimLength, vector::zero)
);
gradPexact.replace(2, P/readScalar(tracerDict.lookup("scaleHeight")));
gradPexact.write();
*/
