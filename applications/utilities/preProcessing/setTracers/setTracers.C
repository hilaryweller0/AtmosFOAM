#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "advectable.H"
#include "tracerField.H"
#include "velocityField.H"
#include "zeroVelocityField.H"
#include "linear.H"
using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    Foam::argList::addOption
    (
        "name",
        "tracerName",
        "specify a non-default tracer name"
    );
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
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

    // list of the names of the tracers
    word tracerName = args.optionFound("name") ?
                      args.optionRead<word>("name") :
                      "T";
    const wordList Tnames
    (
        mesh.solution().lookupOrDefault<wordList>
        (
            "tracerNames",
            wordList(1, tracerName)
        )
    );

    for(label iT = 0; iT < Tnames.size(); iT++)
    {
        tracerName = Tnames[iT];

        Info << "Reading " << tracerName << " if it exists" << endl;
        volScalarField T
        (
            IOobject
            (
                tracerName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            mesh,
            dimensionedScalar(tracerName, dimless, scalar(0))
        );

        const word tracerDictName = tracerName+"InitDict";
        Info<< "Reading initial conditions from " << tracerDictName << endl;

        IOdictionary tracerDict
        (
            IOobject
            (
                tracerDictName,
                runTime.constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        IOdictionary vDict(IOobject("vDict", runTime.constant(), mesh));

        Info << "Setting tracer" << endl;
        autoPtr<velocityField> velocityField(zeroVelocityField::New(vDict));
        autoPtr<tracerField> tracer
        (
            tracerField::New
            (
                tracerDict,
                dynamic_cast<advectable&>(velocityField())
            )
        );
        tracer->applyTo(T);
        T.write();
        Info << "writing " << tracerName << "f for time " << runTime.timeName()
             << endl;
    }

    return EXIT_SUCCESS;
}


