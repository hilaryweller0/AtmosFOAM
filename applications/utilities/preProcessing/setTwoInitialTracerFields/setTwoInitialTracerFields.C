#include "fvCFD.H"
#include "noAdvection.H"
#include "tracerField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading r_init" << endl;
    volScalarField r_init
    (
        IOobject("r_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info << "Reading rho_init" << endl;
    volScalarField rho_init
    (
        IOobject("rho_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info << "Reading or creating tracer field rhof_init" << endl;
    surfaceScalarField rhof_init
    (
        IOobject("rhof_init", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(rho_init)
    );
    
    Info << "Creating r" << endl;
    volScalarField r
    (
        IOobject("r", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating rho" << endl;
    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );

    Info << "Creating rho1" << endl;
    volScalarField rho1
    (
        IOobject("rho1", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );

    Info << "Creating rho2" << endl;
    volScalarField rho2
    (
        IOobject("rho2", runTime.timeName(), mesh, IOobject::NO_READ),
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
    
    IOdictionary rDict
    (
        IOobject
        (
            "rDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> tracer(tracerField::New(tracerDict, velocityField));
    
    autoPtr<tracerField> rVal(tracerField::New(rDict, velocityField));
    
    Info << "writing rho for time " << runTime.timeName() << endl;
    tracer->applyTo(rho);
    rho.write();

    Info << "writing rho1 for time " << runTime.timeName() << endl;
    tracer->applyTo(rho1);
    rho1.write();

    Info << "writing rho2 for time " << runTime.timeName() << endl;
    rho2.write();

    Info << "writing rhof for time " << runTime.timeName() << endl;
    tracer->applyTo(rhof);
    rhof.write();
    
    Info << "writing r" << endl;
    rVal->applyTo(r);
    r.write();

    return EXIT_SUCCESS;
}


