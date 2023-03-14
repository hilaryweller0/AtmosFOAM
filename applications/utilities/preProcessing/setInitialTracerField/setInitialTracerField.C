#include "fvCFD.H"
#include "noAdvection.H"
#include "tracerField.H"

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

    const word tracerName = args.optionFound("name") ?
                              args.optionRead<word>("name") :
                              "T_analytic";

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
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> tracer(tracerField::New(tracerDict, velocityField));

    Info << "setting " << tracerName << " for time " << runTime.timeName()
         << endl;
    tracer->applyTo(T);
    Info << "writing  " << tracerName << endl;
    T.write();

    return EXIT_SUCCESS;
}


