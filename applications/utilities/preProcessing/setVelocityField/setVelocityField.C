#include "fvCFD.H"
#include "velocityField.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
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
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<velocityField> v(velocityField::New(dict));

    forAll(timeDirs, timeI)
    {
        Info << "writing phi for time " << runTime.timeName() << endl;

        v->applyTo(phi);
        phi.write();

        U = fvc::reconstruct(phi);
        Uf = linearInterpolate(U);
        Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        U.write();
        Uf.write();
        
        // Write out the velocity magnitude on faces
        surfaceScalarField magU("magU", mag(Uf));
        magU.write();
        // Write out the face areas
        surfaceScalarField magSf("magSf", mesh.magSf());
        magSf.write();
    }

    return EXIT_SUCCESS;
}


