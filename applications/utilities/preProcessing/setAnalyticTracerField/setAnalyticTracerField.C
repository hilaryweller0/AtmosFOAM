#include "fvCFD.H"
#include "advectable.H"
#include "tracerField.H"
#include "velocityField.H"

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "tracerDict",
        "dictName",
        "specify non-default dictionary name for the tracer (in system)"
    );
    Foam::argList::addOption
    (
        "velocityDict",
        "dictName",
        "specify non-default dictionary name for the velocity (in system)"
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
    timeSelector::addOptions();
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
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

    const word tracerName = args.optionFound("name") ?
                              args.optionRead<word>("name") :
                              "T_analytic";

    Info << "Reading " << tracerName << "_init" << endl;
    volScalarField T_init
    (
        IOobject
        (
            tracerName+"_init",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    Info << "Creating tracer field " << tracerName << endl;
    volScalarField T
    (
        IOobject(tracerName, runTime.timeName(), mesh, IOobject::NO_READ),
        T_init
    );

    const word tracerDictName = args.optionFound("tracerDict") ?
                                args.optionRead<word>("tracerDict") :
                                "tracerFieldDict";
    Info<< "Reading initial conditions from" << tracerDictName << endl;

    IOdictionary tracerDict
    (
        IOobject
        (
            tracerDictName,
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const word velocityDictName = args.optionFound("velocityDict") ?
                                  args.optionRead<word>("velocityDict") :
                                  "velocityFieldDict";
    Info<< "Reading initial conditions from" << velocityDictName << endl;
    IOdictionary velocityDict
    (
        IOobject
        (
            velocityDictName,
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<velocityField> velocityField(velocityField::New(velocityDict));
    autoPtr<tracerField> tracer
    (
        tracerField::New
        (
            tracerDict,
            dynamic_cast<advectable&>(velocityField())
        )
    );

    Info << "Reading or creating tracer field Tf_init" << endl;
    surfaceScalarField Tf_init
    (
        IOobject
        (
            tracerName+"f_init",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        linearInterpolate(T_init)
    );

    Info << "Creating Tf_analytic" << endl;
    surfaceScalarField Tf
    (
        IOobject(tracerName+"f_analytic", runTime.timeName(), mesh),
        Tf_init
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "writing " << tracerName << " for time " << runTime.timeName()
             << endl;
        tracer->applyTo(T);
        T.write();

        Info << "writing " << tracerName << "f_analytic for time "
             << runTime.timeName() << endl;
        tracer->applyTo(Tf);
        Tf.write();
    }

    return EXIT_SUCCESS;
}


