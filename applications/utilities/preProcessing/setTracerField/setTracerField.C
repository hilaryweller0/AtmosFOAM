#include "fvCFD.H"
#include "tracerField.H"
#include "velocityField.H"

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

    IOdictionary velocityDict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<tracerField> tracer(tracerField::New(tracerDict));
//    autoPtr<velocityField> velocity(velocityField::New(velocityDict));

    do
    {
        Info << "writing T for time " << runTime.timeName() << endl;

        tracer->applyTo(T);
        T.write();
    }
    while (runTime.loop());

    return EXIT_SUCCESS;
}


