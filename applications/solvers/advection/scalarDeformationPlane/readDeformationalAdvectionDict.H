Info<< "Reading deformationalAdvectionDict\n" << endl;

IOdictionary initDict
(
    IOobject
    (
        "deformationalAdvectionDict",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Constants in the specification of the velocity field
const scalar endTime(dimensionedScalar(initDict.lookup("endTime")).value());
const scalar psiScale(readScalar(initDict.lookup("psiScale")));
const vector L(dimensionedVector(initDict.lookup("L")).value());

