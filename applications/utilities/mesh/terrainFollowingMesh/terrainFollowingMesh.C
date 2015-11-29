#include "fvCFD.H"
#include "Mountain.H"
#include "BTF.H"

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

    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

    autoPtr<Mountain> mountain(Mountain::New(dict));
    BTF btf(mountain, dict);
    forAll(newPoints, pointIdx)
    {
        newPoints[pointIdx] = btf.computationalToPhysical(newPoints[pointIdx]);
    }

    newPoints.write();
}
