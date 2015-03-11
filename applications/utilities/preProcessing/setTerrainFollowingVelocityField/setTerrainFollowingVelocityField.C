#include "fvCFD.H"
#include "Mountain.H"
#include "VelocityProfile.H"
#include "SchaerCosVelocityProfile.H"
#include "VelocityField.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    Info << "Reading velocityFieldDict" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar u0(readScalar(dict.lookup("maxVelocity")));
    const scalar H(readScalar(dict.lookup("domainHeight")));

    const SchaerCosMountain mountain(dict);
    const SchaerCosVelocityProfile profile(mountain, u0, H);
    const VelocityField velocityField(profile);

    Info << "Creating velocity field Uf" << endl;
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    velocityField.applyTo(Uf);

    Uf.write();

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    Info << "Creating velocity field U" << endl;
    velocityField.applyTo(U);

    U.write();
}
