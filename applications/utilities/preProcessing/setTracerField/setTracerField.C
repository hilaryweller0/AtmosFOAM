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
        "tracerDict",
        "dictName",
        "specify non-default dictionary name for the tracer (in constant)"
    );
    Foam::argList::addOption
    (
        "velocityDict",
        "dictName",
        "specify non-default dictionary name for the velocity (in constant)"
    );
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
    Foam::argList::addBoolOption
    (
        "increment",
        "add the prescribed field as an increment to the existing field"
    );
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.name() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.name(), runTime, IOobject::MUST_READ
        )
    );

    const word tracerName = args.optionFound("name") ?
                              args.optionRead<word>("name") :
                              "T_analytic";

    Info << "Reading " << tracerName << " if it exists" << endl;
    volScalarField T
    (
        IOobject
        (
            tracerName,
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar(tracerName, dimless, scalar(0))
    );

    const word tracerDictName = args.optionFound("tracerDict") ?
                                args.optionRead<word>("tracerDict") :
                                "tracerFieldDict";
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

    const word velocityDictName = args.optionFound("velocityDict") ?
                                  args.optionRead<word>("velocityDict") :
                                  "velocityFieldDict";
    Info<< "Reading initial conditions from " << velocityDictName << endl;
    IOdictionary velocityDict
    (
        IOobject
        (
            velocityDictName,
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    Info << "Setting velocityField" << endl;
    autoPtr<velocityField> velocityField(velocityField::New(velocityDict));
    if (dynamic_cast<advectable*>(&velocityField()) == nullptr)
    {
        WarningIn("setTracerField")
            << "Cannot define analytic solution for tracer for non-advectable "
            << "velocityField. Using zero velocity field" << endl;
    
        velocityField = new zeroVelocityField(velocityDict);
    }
    Info << "Setting tracer" << endl;
    autoPtr<tracerField> tracer
    (
        tracerField::New
        (
            tracerDict,
            dynamic_cast<advectable&>(velocityField())
        )
    );

    Info << "Reading or creating tracer field " << tracerName << "f" << endl;
    surfaceScalarField Tf
    (
        IOobject
        (
            tracerName+"f",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        linearInterpolate(T)
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "writing " << tracerName << " for time " << runTime.name()
             << endl;
        if (args.optionFound("increment"))
        {
            T.oldTime();
            tracer->applyTo(T);
            T == T + T.oldTime();
        }
        else
        {
            tracer->applyTo(T);
        }
        T.write();

        Info << "writing " << tracerName << "f for time "
             << runTime.name() << endl;
        tracer->applyTo(Tf);
        Tf.write();
    }

    return EXIT_SUCCESS;
}


