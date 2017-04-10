#include "fvCFD.H"
#include "noAdvection.H"
#include "tracerField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading rho_init" << endl;
    volScalarField rho_init
    (
        IOobject("rho_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );
    
    Info << "Reading q_init" << endl;
    volScalarField q_init
    (
        IOobject("q_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info << "Reading or creating tracer field rhof_init" << endl;
    surfaceScalarField rhof_init
    (
        IOobject("rhof_init", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(rho_init)
    );

    Info << "Creating q" << endl;
    volScalarField q
    (
        IOobject("q", runTime.timeName(), mesh, IOobject::NO_READ),
        q_init
    );

    Info << "Creating rho" << endl;
    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );

    Info << "Creating rhof" << endl;
    surfaceScalarField rhof
    (
        IOobject("rhof", runTime.timeName(), mesh, IOobject::NO_READ),
        rhof_init
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
    
    IOdictionary qDict
    (
        IOobject
        (
            "qFieldDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> tracer(tracerField::New(tracerDict, velocityField));
    autoPtr<tracerField> qtracer(tracerField::New(qDict, velocityField));

    Info << "writing q for time " << runTime.timeName() << endl;
    qtracer->applyTo(q);
    q.write();

    Info << "writing rho for time " << runTime.timeName() << endl;
    tracer->applyTo(rho);
    rho.write();

    Info << "writing rhof for time " << runTime.timeName() << endl;
    tracer->applyTo(rhof);
    rhof.write();

    return EXIT_SUCCESS;
}


