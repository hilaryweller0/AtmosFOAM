// Create the mesh displacement potential, Phi
volScalarField Phi
(
    IOobject("Phi",runTime.timeName(),mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
Phi.oldTime();

volVectorField gradPhi
(
    IOobject("gradPhi", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::grad(Phi),
    "calculated"
);

// determinant of the Hessian of Phi
volScalarField detHess
(
    IOobject("detHess", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("", dimless, scalar(1))
);
volTensorField Hessian("Hessian", fvc::grad(gradPhi));
volScalarField del2Phi("del2Phi", fvc::laplacian(Phi));
    
volScalarField monitor
(
    IOobject("monitor", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

volScalarField equidist
(
    IOobject("equidist", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    monitor*detHess
);

volScalarField monitorNew("monitorNew", monitor);

volScalarField monitorR
(
    IOobject("monitorR", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    rMesh,
    monitor.dimensions(),
    monitor.internalField(),
    monitor.boundaryField()
);

