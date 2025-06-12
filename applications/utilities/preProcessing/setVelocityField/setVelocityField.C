#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "velocityField.H"
#include "fvcReconstruct.H"
#include "linear.H"
using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "dict", "dictName", "specify the dictionary name (in system)"
    );
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.name() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.name(), runTime, IOobject::MUST_READ
        )
    );

    Info << "Reading in density, rho, if present" << endl;
    volScalarField rho
    (
        IOobject("rho", runTime.name(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimless, scalar(1))
    );
    surfaceScalarField rhof
    (
        IOobject("rhof", runTime.name(), mesh, IOobject::READ_IF_PRESENT),
        linearInterpolate(rho)
    );

    Info << "Creating flux field phi, phiv, U and Uf" << endl;
    surfaceScalarField phi
    (
        IOobject("phi", runTime.name(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimVelocity*dimArea*rho.dimensions(), scalar(0))
    );
    surfaceScalarField phiv
    (
        IOobject("phiv", runTime.name(), mesh, IOobject::READ_IF_PRESENT),
        phi/rhof
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::reconstruct(phiv)
    );

    // Read Uf if present, otherwise create and write
    surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U)
    );

    const word dictName = args.optionFound("dict") ?
                          args.optionRead<word>("dict") :
                          "velocityFieldDict";
    Info<< "Reading initial conditions from " << dictName << endl;
    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<velocityField> v(velocityField::New(dict));

    forAll(timeDirs, timeI)
    {
        Info << "writing phi for time " << runTime.name() << endl;

        v->applyTo(phiv);
        phi = phiv*rhof;
        phiv.write();
        phi.write();

        U = fvc::reconstruct(phiv);
        Uf = linearInterpolate(U);
        Uf += (phiv - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        Info << "writing U and Uf for time " << runTime.name() << endl;
        U.write();
        Uf.write();
    }

    return EXIT_SUCCESS;
}


