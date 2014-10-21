#include "fvCFD.H"

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("defineAtFaces", "Define the velocity vector at cell faces instead of cell centres");
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

    if (args.options().found("defineAtFaces")) {
        surfaceVectorField Uf
        (
            IOobject("Uf", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
            "fixedValue"
        );

        Info << "Creating velocity field U" << endl;
        forAll(Uf, cellI) {
            const point& face = mesh.Cf()[cellI];

            if (face.z() > z1 && face.z() < z2) {
                Uf[cellI] = vector(u0*pow((Foam::sin(M_PI/2*(face.z()-z1)/(z2-z1))),2), 0, 0);
            } else if (face.z() >= z2) {
                Uf[cellI] = vector(u0, 0, 0);
            }
        }

        // TODO: define Uf at inlet and outlet boundaries

   	    Uf.write();
    } else {
        // TODO get some reuse between these two blocks!
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
}
