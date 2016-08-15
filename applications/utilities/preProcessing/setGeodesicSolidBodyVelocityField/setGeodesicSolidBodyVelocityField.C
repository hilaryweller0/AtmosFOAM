#include "fvCFD.H"
#include "polarPoint.H"
#include "sphericalVector.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    dimensionedScalar radius("radius", dimLength, 6.3712e6);
    dimensionedScalar timeScale = runTime.endTime();

    forAll(Uf, faceI)
    {
        const polarPoint& polarp = convertToPolar(mesh.Cf()[faceI]);
        const scalar lat = polarp.lat();
        const scalar lon = polarp.lon();

        scalar alpha0 = 0;
        dimensionedScalar u0 = 2 * M_PI * radius / timeScale;

        dimensionedScalar u = u0 * (Foam::cos(alpha0) * Foam::cos(lat) + Foam::sin(alpha0) * Foam::cos(lon) * Foam::sin(lat));
        dimensionedScalar v = -u0 * Foam::sin(alpha0) * Foam::sin(lon);

        sphericalVector localWind(u.value(), v.value(), 0);
        sphericalVector p(mesh.Cf()[faceI]);

        Uf[faceI] = localWind.toCartesian(p);
    }
    Uf.write();

    return EXIT_SUCCESS;
}


