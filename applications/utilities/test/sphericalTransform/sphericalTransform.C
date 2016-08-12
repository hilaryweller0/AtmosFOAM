#include "fvCFD.H"

int main(int argc, char *argv[])
{
    /*const surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );*/

    vector x(0,0,1);
    const vector kHat(0,0,1);
    vector rHat = x / mag(x);
    if (mag(x) < VSMALL) {
        Info << "Centre of the Earth?!" << endl;
        return EXIT_FAILURE;
    }
    if (mag(kHat ^ rHat) < VSMALL)
    {
        if ((kHat & rHat) > VSMALL)
        {
            Info << "North Pole!" << endl;
            return EXIT_SUCCESS;
        }
        else
        {
            Info << "South Pole!" << endl;
            return EXIT_SUCCESS;
        }
    }
    vector latHat = rHat ^ (kHat ^ rHat);
    vector lonHat = latHat ^ rHat;

    tensor T(latHat, lonHat, rHat);

    vector uLocal(3, 6, 1);
    vector uGlobal = T.inv() & uLocal;

    Info << "x " << x << endl;
    Info << "lon/lat/r " << lonHat << " " << latHat << " " << rHat << endl;
    Info << "uLocal " << uLocal << endl;
    Info << "uGlobal " << uGlobal << endl;

    return EXIT_SUCCESS;
}
