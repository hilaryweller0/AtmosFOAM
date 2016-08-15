#include "fvCFD.H"
#include "polarPoint.H"
#include "sphericalVector.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    IOdictionary velocityFieldDict
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

    const dimensionedScalar radius("radius", dimLength, velocityFieldDict.lookupOrDefault<scalar>("radius", scalar(6.3712e6)));
    const scalar alpha0(velocityFieldDict.lookupOrDefault<scalar>("tilt", scalar(0)));
    const dimensionedScalar timeScale = runTime.endTime();

    forAll(Uf, faceI)
    {
        const polarPoint& polarp = convertToPolar(mesh.Cf()[faceI]);
        const scalar lat = polarp.lat();
        const scalar lon = polarp.lon();

        dimensionedScalar u0 = 2 * M_PI * radius / timeScale;

        dimensionedScalar u = u0 * (Foam::cos(alpha0) * Foam::cos(lat) + Foam::sin(alpha0) * Foam::cos(lon) * Foam::sin(lat));
        dimensionedScalar v = -u0 * Foam::sin(alpha0) * Foam::sin(lon);

        sphericalVector localWind(u.value(), v.value(), 0);
        sphericalVector p(mesh.Cf()[faceI]);

        Uf[faceI] = localWind.toCartesian(p);
    }
    Uf.write();

    Info << "Calculating the divergence field to check that it is zero" << endl;
    volScalarField divu("divu", fvc::div(Uf & mesh.Sf()));
    divu.write();

    Info << "Calculating the radial wind component w to check that it is zero" << endl;
    surfaceScalarField w("w", Uf & (mesh.Cf()/mag(mesh.Cf())));
    w.write();

    return EXIT_SUCCESS;
}


