#include "fvCFD.H"
#include "cellSet.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    forAll(mesh.cells(), cellI)
    {
        const point& p = mesh.C()[cellI];
        forAll(mesh.cellCells()[cellI], neiCellI)
        {
            const point& o = mesh.C()[mesh.cellCells()[cellI][neiCellI]];
            Info << mag(p - o) << endl;
        }
    }

    return EXIT_SUCCESS;
}
