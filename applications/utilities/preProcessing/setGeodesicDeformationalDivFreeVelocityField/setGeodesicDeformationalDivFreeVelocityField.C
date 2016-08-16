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

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", cmptMultiply(dimVelocity, dimArea), vector(0,0,0)),
       "fixedValue"
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
        Uf = phi * mesh.Sf() / mag(mesh.Sf());
        Uf.write();

        runTime.loop();
    }

    return EXIT_SUCCESS;
}


