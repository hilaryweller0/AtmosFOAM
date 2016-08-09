#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating initial tracer field T" << endl;
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("T", dimless, scalar(0)),
        "zeroGradient"
    );

    // section 2.2.1 doi:10.5194/gmd-5-887-2012
    dimensionedScalar radius("radius", dimLength, 6.3712e6);
    dimensionedScalar hmax("hmax", dimless, 0.95);
    dimensionedScalar b("b", dimless, 5);

    scalar lon1 = 5.0*M_PI/6.0;
    scalar lon2 = 7.0*M_PI/6.0;
    scalar lat = 0;

    const dimensionedVector centre1("centre1", dimLength, point(
            radius.value() * Foam::cos(lat) * Foam::cos(lon1),
            radius.value() * Foam::cos(lat) * Foam::sin(lon1),
            radius.value() * Foam::sin(lat)
    ));

    const dimensionedVector centre2("centre2", dimLength, point(
            radius.value() * Foam::cos(lat) * Foam::cos(lon2),
            radius.value() * Foam::cos(lat) * Foam::sin(lon2),
            radius.value() * Foam::sin(lat)
    ));

    T = hmax * exp(-b*magSqr(mesh.C() - centre1)/sqr(radius)) +
        hmax * exp(-b*magSqr(mesh.C() - centre2)/sqr(radius));

    T.write();

    return EXIT_SUCCESS;
}


