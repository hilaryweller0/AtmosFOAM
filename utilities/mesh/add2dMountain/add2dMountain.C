/*---------------------------------------------------------------------------*\
Date started: 10/05/2013

Moves the grid points.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "mountainTypes.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary initDict
    (
        IOobject
        (
            "add2dMountainDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // New point locations layered over the mountain
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );
    
    // Mountain function type
    const word mountainName(initDict.lookup("mountainType"));
    // Declare function pointers to smooth and hi res mountains
    scalar (*smoothMountain)(scalar, scalar, scalar) = 
        mountainName == "ScharExp" ? &ScharExp :
        mountainName == "ScharCos" ? &ScharCosSmooth :
        mountainName == "AgnessiWitch" ? &AgnessiWitch : NULL;
        
    scalar (*fineMountain)(scalar, scalar) = 
        mountainName == "ScharExp" ? &ScharCos :
        mountainName == "ScharCos" ? &ScharCos :
        mountainName == "AgnessiWitch" ? &flatMountain : NULL;

    if 
    (
        mountainName != "ScharExp" && mountainName != "ScharCos"
     && mountainName != "AgnessiWitch"
    )
    {
        FatalErrorIn("ScharMountain")
      <<"mountainType should be one of ScharExp, ScharCos or AgnessiWitch"
         << " not " << mountainName << exit(FatalError);
    }
    
    // Get which coord system to use
    enum coordSysType{BTF, HTF, SLEVE, NONE};
    const word coordSysName(initDict.lookup("coordSys"));
    const coordSysType coordSys = coordSysName == "BTF" ? BTF :
                                  coordSysName == "HTF" ? HTF :
                                  coordSysName == "SLEVE" ? SLEVE : NONE;
    if (coordSys == NONE)
    {
        FatalErrorIn("ScharMountain")
            << "coordSys must be one of BTF, HTF or SLEVE. Not "
            << coordSysName << exit(FatalError);
    }
    
    // Declare and read in constants
    // top of movement of levels
    const scalar zt(readScalar(initDict.lookup("zt")));
    // horizontal mountain scale
    const scalar a(readScalar(initDict.lookup("a")));
    // Maximum mountain height
    const scalar hm(readScalar(initDict.lookup("hm")));

    // horizontal scale in Schar mountain width
    const scalar lambda(initDict.lookupOrDefault<scalar>("lambda", scalar(0)));
    if
    (
        !initDict.found("lambda")
     && (mountainName == "ScharExp" || mountainName == "ScharCos")
    )
    {
        FatalErrorIn("ScharMountain")
            << "if mountain type is Schar, must specify lambda"
            << exit(FatalError);
    }
    
    // Calculate new points for BTF
    switch (coordSys)
    {
    case BTF:
        forAll(newPoints, ip)
        {
            scalar x = newPoints[ip].x();
            scalar z = newPoints[ip].z();
            if (z < zt)
            {
                scalar h = smoothMountain(x,a,hm) * fineMountain(x,lambda);
                newPoints[ip].z() += h*(1 - z/zt);
            }
        }
        break;
        
    case HTF:
    {
        // scale height for HTF
        const scalar zh(readScalar(initDict.lookup("HTF_scaleHeight")));
        forAll(newPoints, ip)
        {
            scalar x = newPoints[ip].x();
            scalar z = newPoints[ip].z();
            scalar h = smoothMountain(x,a,hm) * fineMountain(x,lambda);
            if (z < zh && z < zt)
            {
                newPoints[ip].z() += h*pow(Foam::cos(0.5*M_PI*z/zh),6);
            }
        }
    }break;
    
    case SLEVE:
    {
        // scale heights for coarse and fine mountains and the exponent
        FixedList<scalar,2> s, h, b;
        s[0] = readScalar(initDict.lookup("SLEVE_scaleCoarse"));
        s[1] = readScalar(initDict.lookup("SLEVE_scaleFine"));
        const scalar n(readScalar(initDict.lookup("SLEVE_exponent")));
        forAll(newPoints, ip)
        {
            scalar z = newPoints[ip].z();
            if (z < zt)
            {
                scalar x = newPoints[ip].x();
                h[0] = 0.5*smoothMountain(x,a,hm);
                h[1] = 2*h[0]*fineMountain(x,lambda) - h[0];
                for(int i = 0; i < 2; i++)
                {
                    b[i] = Foam::sinh
                           (
                               Foam::pow(zt/s[i], n)
                             - Foam::pow(z/s[i], n)
                           )
                           /Foam::sinh(Foam::pow(zt/s[i],n));
                    newPoints[ip].z() += h[i]*b[i];
                }
            }
        }
    }break;
    
    default:
        FatalErrorIn("ScharMountain")
            << "coordSys must be one of BTF, HTF or SLEVE. Not "
            << coordSysName << exit(FatalError);
    }

    newPoints.write();
}
