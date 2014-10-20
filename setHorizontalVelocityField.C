#include "fvCFD.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    Info << "Reading velocityFieldDict" << endl;

    IOdictionary initDict
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

    const scalar u0(readScalar(initDict.lookup("maxVelocity")));
    const scalar z1(readScalar(initDict.lookup("maxVelocityHeight")));
    const scalar z2(readScalar(initDict.lookup("zeroVelocityHeight")));

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0,0,0)),
        "zeroGradient"
    );

    Info << "Creating velocity field U" << endl;
    forAll(U, cellI) {
        const point& c = mesh.C()[cellI];

        if (c.z() > z1 && c.z() < z2) {
            U[cellI] = vector(u0*pow((Foam::sin(M_PI/2*(c.z()-z1)/(z2-z1))),2), 0, 0);
        } else if (c.z() >= z2) {
            U[cellI] = vector(u0, 0, 0);
        }
    }

   	U.write();
}
