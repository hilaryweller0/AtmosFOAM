#include "Time.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "mountain.H"

using namespace Foam;

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
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar tol(readScalar(dict.lookup("tolerance")));
    
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

    autoPtr<mountain> m(mountain::New(dict));

	forAll(newPoints, pointIdx)
	{
        scalar h = m->heightAt(newPoints[pointIdx]).value();
        if (newPoints[pointIdx].z() < h + tol)
        {
            newPoints[pointIdx].z() = h;
        }
    }

    newPoints.write();
}
