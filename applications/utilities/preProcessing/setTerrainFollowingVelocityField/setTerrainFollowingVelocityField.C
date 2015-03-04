#include "fvCFD.H"
#include "Mountain.H"

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    Info << "Creating velocity field Uf" << endl;
    // TODO

    SchaerCosMountain mountain(1.0, 2.0, 3.0);
    
    /*
     * mountain (to get height for a given x)
     * u0 -- base value of u
     */

    Uf.write();
}
