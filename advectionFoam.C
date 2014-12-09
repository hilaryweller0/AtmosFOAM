#include "fvCFD.H"

volScalarField readT(Time& runTime, fvMesh& mesh)
{
    Info<< "Reading field T\n" << endl;

    return volScalarField
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}

volScalarField calculateDivUf(Time& runTime, fvMesh& mesh, surfaceScalarField& phi)
{
    Info<< "Calculating divergence of phi\n" << endl;

    return volScalarField
    (
        IOobject(
            "divUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi)
    );
}

surfaceScalarField readOrCalculatePhi(argList& args, Time& runTime, fvMesh& mesh)
{
    if (args.options().found("usePhi"))
    {
        Info<< "Reading field phi\n" << endl;

        return surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Reading field Uf\n" << endl;

        surfaceVectorField Uf
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Calculating face flux field phi\n" << endl;

        return surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Uf & mesh.Sf()
        );
    }
}

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("usePhi", "use 0/phi rather than calculating it from 0/Uf");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField T = readT(runTime, mesh);
    surfaceScalarField phi = readOrCalculatePhi(args, runTime, mesh);
    volScalarField divUf = calculateDivUf(runTime, mesh, phi);
    divUf.write();

    Info<< "\nCalculating advection\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int corr=0; corr < 3; corr++)
        {
            T = T.oldTime() - runTime.deltaT() * 0.5 *
            (
                fvc::div(phi, T) + fvc::div(phi, T.oldTime())
            );
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
