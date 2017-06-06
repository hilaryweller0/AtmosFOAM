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
    
    Info << "Reading P_init" << endl;
    volScalarField P_init
    (
        IOobject("P_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

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
    
    Info << "Creating T" << endl;
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ),
        T_init
    );
    
    Info << "Creating P" << endl;
    volScalarField P
    (
        IOobject("P", runTime.timeName(), mesh, IOobject::NO_READ),
        P_init
    );
    
    Info << "Creating rho" << endl;
    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );
    
    Info << "Creating rho_a" << endl;
    volScalarField rho_a
    (
        IOobject("rho_a", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );
    
    Info << "Creating rho_b" << endl;
    volScalarField rho_b
    (
        IOobject("rho_b", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );

    Info << "Creating q1" << endl;
    volScalarField q1
    (
        IOobject("q1", runTime.timeName(), mesh, IOobject::NO_READ),
        q_init
    );

    Info << "Creating q2" << endl;
    volScalarField q2
    (
        IOobject("q2", runTime.timeName(), mesh, IOobject::NO_READ),
        rho_init
    );
    
    Info << "Creating S" << endl;
    volScalarField S
    (
        IOobject("S", runTime.timeName(), mesh, IOobject::NO_READ),
        q_init
    );

    Info << "Creating rhof" << endl;
    surfaceScalarField rhof
    (
        IOobject("rhof", runTime.timeName(), mesh, IOobject::NO_READ),
        rhof_init
    );

    IOdictionary rhoDict
    (
        IOobject
        (
            "rhoDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary rhoaDict
    (
        IOobject
        (
            "rhoaDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary qFieldDict
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
    
    IOdictionary tempDict
    (
        IOobject
        (
            "tempDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary PDict
    (
        IOobject
        (
            "PDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> rhoVal(tracerField::New(rhoDict, velocityField));
    
    autoPtr<tracerField> rhoaVal(tracerField::New(rhoaDict, velocityField));
    
    autoPtr<tracerField> qVal(tracerField::New(qFieldDict, velocityField));
    
    autoPtr<tracerField> tempVal(tracerField::New(tempDict, velocityField));
    
    autoPtr<tracerField> PVal(tracerField::New(PDict, velocityField));
    
    Info << "writing rho for time " << runTime.timeName() << endl;
    rhoVal->applyTo(rho);
    rho.write();
    
    Info << "writing rho for time " << runTime.timeName() << endl;
    rhoaVal->applyTo(rho_a);
    rho_a.write();

    Info << "writing rho for time " << runTime.timeName() << endl;
    rhoVal->applyTo(rho_b);
    rho_b.write();

    Info << "writing q1 for time " << runTime.timeName() << endl;
    q1.write();

    Info << "writing q2 for time " << runTime.timeName() << endl;
    qVal->applyTo(q2);
    q2.write();

    Info << "writing S for time " << runTime.timeName() << endl;
    S.write();

    Info << "writing qf for time " << runTime.timeName() << endl;
    rhoVal->applyTo(rhof);
    rhof.write();
    
    Info << "writing T" << endl;
    tempVal->applyTo(T);
    T.write();
    
    Info << "writing P" << endl;
    PVal->applyTo(P);
    P.write();

    return EXIT_SUCCESS;
}


