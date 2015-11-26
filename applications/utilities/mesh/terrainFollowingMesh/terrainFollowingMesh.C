#include "fvCFD.H"
#include "Mountain.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "mountainDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << "Hello world" << endl;
    
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

	forAll(newPoints, pointIdx)
	{
        // TODO
    }

    newPoints.write();
}
