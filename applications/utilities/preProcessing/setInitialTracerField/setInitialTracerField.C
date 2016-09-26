#include "fvCFD.H"
#include "noAdvection.H"
#include "tracerField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating tracer field T" << endl;
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("T", dimless, scalar(0)),
       "fixedValue"
    );

    Info << "Creating tracer field Tf" << endl;
    surfaceScalarField Tf
    (
        IOobject("Tf", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Tf", dimless, scalar(0)),
       "fixedValue"
    );

    IOdictionary tracerDict
    (
        IOobject
        (
            "tracerFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> tracer(tracerField::New(tracerDict, velocityField));

    Info << "writing T for time " << runTime.timeName() << endl;
    tracer->applyTo(T);
    T.write();

    Info << "writing Tf for time " << runTime.timeName() << endl;
    tracer->applyTo(Tf);
    Tf.write();

    return EXIT_SUCCESS;
}


