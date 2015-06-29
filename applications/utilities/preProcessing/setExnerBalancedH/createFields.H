// Read in prognostic variables

Info<< "Creating Exner_init\n" << endl;
volScalarField Exner_init
(
    IOobject("Exner_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Setting Exner\n" << endl;
volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
    Exner_init
);

Info<< "Reading thetaf\n" << endl;
const surfaceScalarField thetaf
(
    IOobject("thetaf", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("U", dimensionSet(0,3,-2,0,0), scalar(0))
);

label groundBC = -1;
label topBC = -1;
forAll(mesh.boundaryMesh(), patchi)
{
    if (mesh.boundaryMesh()[patchi].name() == "top")
    {
//        ExnerBCs[patchi] = "fixedValue";
        topBC = patchi;
    }
    else if (mesh.boundaryMesh()[patchi].name() == "ground")
    {
        groundBC = patchi;
    }
}
if (groundBC == -1 || topBC == -1)
{
    FatalErrorIn("setExnerBalancedH")
        << " no boundary called top or ground. The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);
}

//// Set field which is zero on vertical faces
//surfaceScalarField horiz
//(
//    IOobject("horiz", runTime.constant(), mesh),
//    mesh,
//    dimensionedScalar("", dimless, scalar(1))
//);
//forAll(horiz, faceI)
//{
//    if (mag(mesh.Sf()[faceI] & g.value()) > SMALL) horiz[faceI] = 0.1;
//}

