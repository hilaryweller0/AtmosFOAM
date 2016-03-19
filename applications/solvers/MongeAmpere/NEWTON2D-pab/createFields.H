// Create the mesh displacement potential, Phi
volScalarField Phi
(
    IOobject("Phi",runTime.timeName(),mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Create the perturbation from Phi
volScalarField phi("phi", Phi);

// Create phiBarLaplacian
volScalarField phiBarLaplacian("phiBarLaplacian", fvc::laplacian(Phi));

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
    IOobject("gradPhif", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::interpolate(gradPhi),
    "calculated"
);
gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
            *mesh.Sf()/mesh.magSf();

// Hessian of Phi and its determinant
volTensorField Hessian("Hessian", fvc::grad(gradPhif));
volScalarField detHess
(
    IOobject("detHess", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("", dimless, scalar(1))
);

volScalarField detHessR
(
    IOobject("detHessR", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
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

// the volume ratio between old and new meshes (the geometric Hessian)
volScalarField volRatio
(
    IOobject("volRatio", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimless, wordList(5,"zeroGradient")
);
volRatio.internalField() = 1;
volRatio.correctBoundaryConditions();

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
    IOobject("monitorR", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    monitorFunc().map(rMesh, monitor)
);
monitorNew.internalField() = monitorR.internalField();

// The mean equidistribution, c
dimensionedScalar equiDistMean = fvc::domainIntegrate(detHess)
                                 /fvc::domainIntegrate(1/monitorNew);

dimensionedScalar equiDistMeanR = fvc::domainIntegrate(detHessR)
                                 /fvc::domainIntegrate(1/monitorR);
Info << "Dimension of monitor    = " << monitor.dimensions() << endl;
Info << "Dimension of monitorNew = " << monitorNew.dimensions() << endl;
Info << "Dimension of monitorR   = " << monitorR.dimensions() << endl;


equiDistMeanR = equiDistMean;

// The equidistribution
volScalarField equiDist
(
    IOobject("equiDist", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    monitorR*detHess
);

// The global equidistribution
dimensionedScalar PABem = fvc::domainIntegrate(equiDist)/Vtot;
dimensionedScalar PABe = sqrt(fvc::domainIntegrate(sqr(equiDist - PABem)))
                         /(Vtot*PABem);

// constant divided by monitor function:
volScalarField c_m
(
    IOobject("c_m", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    equiDistMean/monitorNew
);

// constant divided by monitor function:
volScalarField c_mR
(
    IOobject("c_mR", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    equiDistMeanR/monitorR
);

// The surface gradient of c_m on on the computational mesh
surfaceScalarField sngradc_m
(
    IOobject("sngradc_m", runTime.timeName(), mesh),
    fvc::snGrad(c_m)
);


// The surface gradient of c_m on the physical mesh
surfaceScalarField sngradc_mR
(
    IOobject("sngradc_mR", runTime.timeName(), rMesh),
    fvc::snGrad(c_mR)
);

// The divergence of sngradc_mR
volScalarField lapc_m
(
    IOobject("lapc_m", runTime.timeName(), mesh),
    fvc::div(mesh.magSf()*sngradc_m)
);

