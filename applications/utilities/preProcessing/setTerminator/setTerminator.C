#include "Time.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "mathematicalConstants.H"
using namespace Foam;
using namespace Foam::constant::mathematical;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading terminator setup from terminatorDict" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "terminatorDict", mesh.time().constant(), mesh, IOobject::MUST_READ
        )
    );
    const scalar latc(readScalar(initDict.lookup("extinctionLatitude"))*pi/180);
    const scalar lonc(readScalar(initDict.lookup("extinctionLongitude"))*pi/180);
    const scalar Xtotal(readScalar(initDict.lookup("Xtotal")));
    const dimensionedScalar k2(initDict.lookup("k2"));

    // Create fields
    volScalarField k1
    (
        IOobject("k1", "constant", runTime),
        mesh,
        dimensionedScalar(dimensionSet(0,0,-1,0,0), scalar(0))
    );
    volScalarField X
    (
        IOobject("X", runTime.timeName(), runTime, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimless, scalar(0))
    );
    volScalarField X2
    (
        IOobject("X2", runTime.timeName(), runTime, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimless, scalar(0))
    );
    volScalarField X1
    (
        IOobject("X1", runTime.timeName(), runTime, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimless, scalar(0))
    );
    volScalarField Xsum
    (
        IOobject("Xsum", runTime.timeName(), runTime, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar(dimless, Xtotal)
    );

    // Calculate the reaction rate k1
    volScalarField x = mesh.C().component(0);
    volScalarField y = mesh.C().component(1);
    volScalarField z = mesh.C().component(2);
    volScalarField R = sqrt(sqr(x) + sqr(y) + sqr(z));
    k1 = max
    (
        k1,
        1/dimensionedScalar(dimTime, scalar(1))*
        (
            Foam::sin(latc)*z/R
          + Foam::cos(latc)*(Foam::cos(lonc)*x/R + Foam::sin(lonc)*y/R)
        )
    );
    k1.write();
    
    // Calculate X, X2 and X1
    volScalarField r = k1/(4*k2);
    volScalarField D = sqrt(sqr(r) + 2*r*Xtotal);
    X = D - r;
    X1 = Xtotal - D + r;
    X2 = 0.5*X1;
    Xsum.write();
    X.write();
    X1.write();
    X2.write();

    return EXIT_SUCCESS;
}


