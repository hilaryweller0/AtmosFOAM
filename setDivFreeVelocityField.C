#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"        

    fvScalarMatrix pEqn
    (
        fvm::laplacian(p) + fvc::div(phi)
    );
    pEqn.solve();

    phi -= pEqn.flux();
    phi.write();
}
