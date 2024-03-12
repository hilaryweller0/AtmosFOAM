#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "GeometricField.H"
#include "processorFvPatchField.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    surfaceVectorField uf
    (
        Foam::IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    unsigned int c = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (mesh.boundaryMesh().types()[patchi] == word("processor"))
        {
            Field<vector>& boundaryField = uf.boundaryFieldRef()[patchi];
            forAll(boundaryField, i)
            {
                if (boundaryField[i][0] < 0)
                {
                    boundaryField[i][0] = -boundaryField[i][0];
                    c++;
                }
            }

            Info << "Corrected " << c << " value(s) in patch " << mesh.boundaryMesh()[patchi].name() << endl;
        }
    }

    uf.write();

    return 0;
}
