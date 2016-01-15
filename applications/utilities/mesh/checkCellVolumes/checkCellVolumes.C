#include "fvCFD.H"
#include "cellSet.H"
#include "Mountain.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    cellSet cells(mesh, "tinyCells", mesh.nCells());
    const scalarField& vols = mesh.cellVolumes();
    scalar maxVolume = GREAT;
    forAll(vols, cellI)
    {
        maxVolume = max(maxVolume, vols[cellI]);
    }
    forAll(vols, cellI)
    {
        if (vols[cellI] < maxVolume/1e12)
        {
            cells.insert(cellI);
        }
    }
    label nCells = returnReduce(cells.size(), sumOp<label>());

    cells.instance() = mesh.pointsInstance();
    cells.write();

    Info << "Found " << nCells << " tiny cells" << endl;
}
