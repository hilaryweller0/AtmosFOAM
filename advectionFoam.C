#include "fvCFD.H"
#include "fvIOoptionList.H"
#include "simpleControl.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    simpleControl simple(mesh);

    Info<< "\nCalculating advection\n" << endl;

    #include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
             ==
                fvOptions(T)
            );
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
