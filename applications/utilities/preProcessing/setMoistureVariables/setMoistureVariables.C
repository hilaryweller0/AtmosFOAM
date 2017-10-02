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

    Info << "Reading rt_init" << endl;
    volScalarField rt_init
    (
        IOobject("rt_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );
    
    Info << "Reading r_init" << endl;
    volScalarField r_init
    (
        IOobject("r_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info << "Reading or creating tracer field rhof_init" << endl;
    surfaceScalarField rhof_init
    (
        IOobject("rhof_init", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(rt_init)
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
    
    Info << "Creating rt" << endl;
    volScalarField rt
    (
        IOobject("rt", runTime.timeName(), mesh, IOobject::NO_READ),
        rt_init
    );

    Info << "Creating rl" << endl;
    volScalarField rl
    (
        IOobject("rl", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );

    Info << "Creating rv" << endl;
    volScalarField rv
    (
        IOobject("rv", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating rl_diag" << endl;
    volScalarField rl_diag
    (
        IOobject("rl_diag", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );

    Info << "Creating rv_diag" << endl;
    volScalarField rv_diag
    (
        IOobject("rv_diag", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating rl_analytic" << endl;
    volScalarField rl_analytic
    (
        IOobject("rl_analytic", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating rv_analytic" << endl;
    volScalarField rv_analytic
    (
        IOobject("rv_analytic", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating rt_analytic" << endl;
    volScalarField rt_analytic
    (
        IOobject("rt_analytic", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );
    
    Info << "Creating S" << endl;
    volScalarField S
    (
        IOobject("S", runTime.timeName(), mesh, IOobject::NO_READ),
        r_init
    );

    Info << "Creating rhof" << endl;
    surfaceScalarField rhof
    (
        IOobject("rhof", runTime.timeName(), mesh, IOobject::NO_READ),
        rhof_init
    );

    IOdictionary rtDict
    (
        IOobject
        (
            "totalMoistureDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary rlDict
    (
        IOobject
        (
            "liquidWaterDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary rvDict
    (
        IOobject
        (
            "waterVapourDict",
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
            "pressureDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const noAdvection velocityField;
    autoPtr<tracerField> rtVal(tracerField::New(rtDict, velocityField));
    
    autoPtr<tracerField> rlVal(tracerField::New(rlDict, velocityField));
    
    autoPtr<tracerField> rvVal(tracerField::New(rvDict, velocityField));
    
    autoPtr<tracerField> tempVal(tracerField::New(tempDict, velocityField));
    
    autoPtr<tracerField> PVal(tracerField::New(PDict, velocityField));
    
    Info << "writing rt for time " << runTime.timeName() << endl;
    rtVal->applyTo(rt);
    rt.write();
    
    Info << "writing rl for time " << runTime.timeName() << endl;
    rlVal->applyTo(rl);
    rl.write();

    Info << "writing rv for time " << runTime.timeName() << endl;
    rvVal->applyTo(rv);
    rv.write();
    
    Info << "writing rl_diag for time " << runTime.timeName() << endl;
    rlVal->applyTo(rl_diag);
    rl_diag.write();

    Info << "writing rv_diag for time " << runTime.timeName() << endl;
    rvVal->applyTo(rv_diag);
    rv_diag.write();

    Info << "writing rl_analytic for time " << runTime.timeName() << endl;
    rlVal->applyTo(rl_analytic);
    rl_analytic.write();
    
    Info << "writing rv_analytic for time " << runTime.timeName() << endl;
    rvVal->applyTo(rv_analytic);
    rv_analytic.write();

    Info << "writing rt_analytic for time " << runTime.timeName() << endl;
    rtVal->applyTo(rt_analytic);
    rt_analytic.write();

    Info << "writing S for time " << runTime.timeName() << endl;
    S.write();

    Info << "writing qf for time " << runTime.timeName() << endl;
    rtVal->applyTo(rhof);
    rhof.write();
    
    Info << "writing T" << endl;
    tempVal->applyTo(T);
    T.write();
    
    Info << "writing P" << endl;
    PVal->applyTo(P);
    P.write();

    return EXIT_SUCCESS;
}


