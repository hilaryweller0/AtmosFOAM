#include "fvCFD.H"
#include "polarPoint.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating wind field Uf" << endl;
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
       "fixedValue"
    );

    // section 2.3 doi:10.5194/gmd-5-887-2012
    dimensionedScalar radius("radius", dimLength, 6.3712e6);
    dimensionedScalar timeScale =runTime.endTime();

    while (runTime.run())
    {
        Info << "writing Uf for time " << runTime.timeName() << endl;
        forAll(Uf, faceI)
        {
            const polarPoint& polarp = convertToPolar(mesh.Cf()[faceI]);
            const scalar lat = polarp.lat();
            const scalar lon = polarp.lon();

            dimensionedScalar t = runTime;
            dimensionedScalar lonPrime = lon - 2 * M_PI * t / timeScale;

            dimensionedScalar u = 10*radius/timeScale * sqr(Foam::sin(lonPrime)) * Foam::sin(2*lat) * Foam::cos(M_PI*t/timeScale) + 2*M_PI*radius/timeScale * Foam::cos(lat);
            dimensionedScalar v = 10*radius/timeScale * Foam::sin(2*lonPrime) * Foam::cos(lat) * Foam::cos(M_PI*t/timeScale);

            scalar latRot = lat+M_PI/2;
            scalar lonRot = lon+M_PI;

            scalar xGlobal = u.value() * Foam::cos(lonRot);
            scalar yGlobal = v.value() * Foam::cos(latRot) + Foam::sin(latRot) * u.value() * Foam::sin(lonRot);
            scalar zGlobal = v.value() * Foam::sin(latRot) - Foam::cos(latRot) * u.value() * Foam::sin(lonRot);

            Uf[faceI] = vector(xGlobal, yGlobal, zGlobal);
        }
        Uf.write();
        runTime.loop();
    }

    return EXIT_SUCCESS;
}


