#include "fvCFD.H"
#include "terrainFollowingTransform.H"

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

    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

    autoPtr<terrainFollowingTransform> transform(terrainFollowingTransform::New(dict));

    forAll(newPoints, pointIdx)
    {
        newPoints[pointIdx] = transform->computationalToPhysical(newPoints[pointIdx]);
    }

    newPoints.write();
}
