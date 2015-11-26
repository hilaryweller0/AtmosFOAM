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
    
    // New point locations layered over the mountain
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

    const word mountainType(dict.lookup("mountainType"));
    Mountain* mountain;
    if (mountainType == "SCHAER_COS")
    {
        mountain = new SchaerCosMountain(dict);
    }
    else if (mountainType == "SCHAER_EXP")
    {
        mountain = new SchaerExpMountain(dict);
    }
    else
    {
        FatalErrorIn("slantMesh")
            << "Unknown mountainType " << mountainType
            << exit(FatalError);
    }

	forAll(newPoints, pointIdx)
	{
        scalar h = mountain->heightAt(newPoints[pointIdx]);
        if (newPoints[pointIdx].z() < h)
        {
            newPoints[pointIdx].z() = h;
        }
    }

    newPoints.write();
}
