// Create the mesh displacement potential, Phi
volScalarField Phi
(
    IOobject("Phi",runTime.timeName(),mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Create the perturbation from Phi
volScalarField phi("phi", Phi);

// The gradient of Phi at cell centres
volVectorField gradPhi
(
    IOobject("gradPhi", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::reconstruct(fvc::snGrad(Phi)*mesh.magSf())
);
gradPhi.boundaryField()
    == (static_cast<volVectorField>(fvc::grad(Phi))).boundaryField();

// The gradient of Phi on faces
surfaceVectorField gradPhif
(
    IOobject("gradPhif", runTime.timeName(), mesh),
    fvc::interpolate(gradPhi),
    "calculated"
);
gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
            *mesh.Sf()/mesh.magSf();

// Hessian of Phi and its determinant
volTensorField Hessian("Hessian", fvc::grad(gradPhif));
volScalarField detHess
(
    IOobject("detHess", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("", dimless, scalar(1))
);

volScalarField detHessR
(
    IOobject("detHessR", runTime.timeName(), rMesh),
    rMesh,
    dimensionedScalar("", dimless, scalar(1))
);


// The reqired matrix A
volTensorField I
(
    IOobject("I", runTime.constant(), mesh),
    mesh,
    dimensionedTensor("I", dimless, tensor::one)
);
volTensorField matrixA("matrixA", -Hessian + I*(1+fvc::laplacian(Phi)));

// The original monitor function on the old mesh
const volScalarField monitor
(
    IOobject("monitor", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("monitor", dimensionSet(0,0,0,0,0), scalar(1)),
    wordList(5,"zeroGradient")
);

// The monitor function on the old mesh calculated from the new mesh
volScalarField monitorNew("monitorNew", monitor);

// The monitor function on the new mesh
volScalarField monitorR
(
    IOobject("monitorR", runTime.timeName(), rMesh),
    monitorFunc().map(rMesh, monitor)
);
monitorNew.internalField() = monitorR.internalField();

// The mean equidistribution, c
dimensionedScalar equiDistMean = fvc::domainIntegrate(detHess)
                                 /fvc::domainIntegrate(1/monitorNew);

// The equidistribution
volScalarField equiDist
(
    IOobject("equiDist", runTime.timeName(), mesh),
    monitorNew*detHess
);

// The global equidistribution
dimensionedScalar PABem = fvc::domainIntegrate(equiDist)/Vtot;
dimensionedScalar PABe = sqrt(fvc::domainIntegrate(sqr(equiDist - PABem)))
                         /(Vtot*PABem);

// constant divided by monitor function:
volScalarField c_m
(
    IOobject("c_m", runTime.timeName(), mesh),
    equiDistMean/monitorNew
);

// constant divided by monitor function on the rMesh
volScalarField c_mR
(
    IOobject("c_mR", runTime.timeName(), rMesh),
    equiDistMean/monitorR
);

// The surface gradient of c_m on the physical mesh
surfaceScalarField sngradc_mR
(
    IOobject("sngradc_mR", runTime.timeName(), rMesh),
    fvc::snGrad(c_mR)
);
// The surface gradient of c_m on on the computational mesh
surfaceScalarField sngradc_m
(
    IOobject("sngradc_m", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("sngradc_m", dimensionSet(0,-1,0,0,0), scalar(0)),
    wordList(5,"calculated")
);
sngradc_m.internalField() = sngradc_mR.internalField();

// The divergence of sngradc_mR
volScalarField lapc_m
(
    IOobject("lapc_m", runTime.timeName(), mesh),
    fvc::div(mesh.magSf()*sngradc_m)
);

