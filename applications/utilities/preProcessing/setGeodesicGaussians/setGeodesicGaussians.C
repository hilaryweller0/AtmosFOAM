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
//    scalar R = 6.3712e6;
    scalar R = 1.0;
    scalar hmax = 0.95;
    scalar b = 5.0;
    scalar lon1 = 5.0*M_PI/6.0;
    scalar lon2 = 7.0*M_PI/6.0;
    scalar lat = 0;

    forAll(T, cellI)
    {
        const point& c = mesh.C()[cellI];

        scalar xi = R * Foam::cos(lat) * Foam::cos(lon1);
        scalar yi = R * Foam::cos(lat) * Foam::sin(lat);
        scalar zi = R * Foam::sin(lat);

        scalar h1 = hmax * Foam::exp(-b * (
                pow(c.x() - xi, 2) +
                pow(c.y() - yi, 2) + 
                pow(c.z() - zi, 2)
        ));

        xi = R * Foam::cos(lat) * Foam::cos(lon2);
        yi = R * Foam::cos(lat) * Foam::sin(lat);
        zi = R * Foam::sin(lat);

        scalar h2 = hmax * Foam::exp(-b * (
                pow(c.x() - xi, 2) +
                pow(c.y() - yi, 2) + 
                pow(c.z() - zi, 2)
        ));
        
        T[cellI] = h1 + h2;
    }

    T.write();

    return EXIT_SUCCESS;
}


