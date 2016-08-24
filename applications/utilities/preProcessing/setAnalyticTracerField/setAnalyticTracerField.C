#include "fvCFD.H"
#include "advectable.H"
#include "tracerField.H"
#include "velocityField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating tracer field T_analytic" << endl;
    volScalarField T
    (
        IOobject("T_analytic", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("T_analytic", dimless, scalar(0)),
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

    autoPtr<velocityField> velocityField(velocityField::New(velocityDict));
    autoPtr<tracerField> tracer
    (
        tracerField::New
        (
            tracerDict,
            dynamic_cast<advectable&>(velocityField())
        )
    );

    do
    {
        Info << "writing T_analytic for time " << runTime.timeName() << endl;
        tracer->applyTo(T);
        T.write();
    } while (runTime.loop());

    return EXIT_SUCCESS;
}


