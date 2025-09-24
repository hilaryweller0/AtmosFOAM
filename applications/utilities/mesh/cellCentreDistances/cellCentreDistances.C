#include "polyMesh.H"
#include "Time.H"
#include "argList.H"
#include "pointField.H"
#include "labelList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    const vectorField& C = mesh.cellCentres();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    
    forAll(mesh.faceNeighbour(), faceI)
    {
        const point& P = C[own[faceI]];
        const point& N = C[nei[faceI]];
        Info << "Distance " << faceI << " " << mag(P - N) << nl;
    }

    return EXIT_SUCCESS;
}
