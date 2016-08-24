#include "fvCFD.H"
#include "velocityField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating flux field phi" << endl;
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("phi", cmptMultiply(dimVelocity, dimArea), scalar(0)),
       "fixedValue"
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::reconstruct(phi)
    );

    // Read Uf if present, otherwise create and write (not used)
    surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(fvc::reconstruct(phi))
    );

    IOdictionary dict
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

    autoPtr<velocityField> v(velocityField::New(dict));

    while (runTime.run())
    {
        Info << "writing phi for time " << runTime.timeName() << endl;

        v->applyTo(phi);
        phi.write();

        U = fvc::reconstruct(phi);
        Uf = linearInterpolate(U);
        Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        U.write();
        Uf.write();

        runTime.loop();
    }

    return EXIT_SUCCESS;
}


