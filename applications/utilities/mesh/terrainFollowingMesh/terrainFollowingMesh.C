#include "fvCFD.H"
#include "terrainFollowingTransform.H"

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
#   include "setRootCase.H"
#   include "createTime.H"
    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
        )
    );

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

    autoPtr<terrainFollowingTransform> transform
    (
        terrainFollowingTransform::New(dict)
    );

    forAll(newPoints, pointIdx)
    {
        newPoints[pointIdx]
             = transform->computationalToPhysical(newPoints[pointIdx]);
    }

    newPoints.write();
}
