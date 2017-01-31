#include "fvCFD.H"
#include "noAdvection.H"
#include "tracerField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading T_init" << endl;
    volScalarField T_init
    (
        IOobject("T_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info << "Reading or creating tracer field Tf_init" << endl;
    surfaceScalarField Tf_init
    (
        IOobject("Tf_init", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(T_init)
    );

    Info << "Creating T1" << endl;
    volScalarField T1
    (
        IOobject("T1", runTime.timeName(), mesh, IOobject::NO_READ),
        T_init
    );

    Info << "Creating T2" << endl;
    volScalarField T2
    (
        IOobject("T2", runTime.timeName(), mesh, IOobject::NO_READ),
        T_init
    );

    Info << "Creating Tf" << endl;
    surfaceScalarField Tf
    (
        IOobject("Tf", runTime.timeName(), mesh, IOobject::NO_READ),
        Tf_init
    );

    IOdictionary tracerDict
    (
        IOobject
        (
            "tracerFieldDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> tracer(tracerField::New(tracerDict, velocityField));

    Info << "writing T1 for time " << runTime.timeName() << endl;
    tracer->applyTo(T1);
    T1.write();

    Info << "writing T2 for time " << runTime.timeName() << endl;
    T2.write();

    Info << "writing Tf for time " << runTime.timeName() << endl;
    tracer->applyTo(Tf);
    Tf.write();

    return EXIT_SUCCESS;
}


