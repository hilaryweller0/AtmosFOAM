#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "mountainTypes.H"
#include "HashTable.H"
#include "fvMeshSubset.H"
#include <cmath>

using namespace Foam::constant::mathematical;

typedef scalar (*SmoothMountainFunction)(scalar, scalar, scalar);
typedef scalar (*FineMountainFunction)(scalar, scalar);

class Mountain
{
    public:
    Mountain(
            const scalar zt,
            const scalar a,
            const scalar hm,
            const scalar lambda,
            const SmoothMountainFunction smoothMountain,
            const FineMountainFunction fineMountain)
        :
            zt(zt),
            a(a),
            hm(hm),
            lambda(lambda),
            smoothMountain(smoothMountain),
            fineMountain(fineMountain)
    {}

    scalar heightAt(const scalar x) const
    {
        return smoothMountain(x,a,hm) * fineMountain(x,lambda);
    }

    private:
    const scalar zt; // top of movement of levels
    const scalar a; // horizontal mountain scale
    const scalar hm; // Maximum mountain height
    const scalar lambda; // horizontal scale in Schar mountain width
    const SmoothMountainFunction smoothMountain;
    const FineMountainFunction fineMountain;
};

int roundUp(int numToRound, int multiple) 
{
   return (numToRound + multiple - 1) / multiple * multiple;
}

void badCoordinateSystem(const word& coordSysName)
{
        FatalErrorIn("add2dMountain")
            << "coordSys must be one of BTF, HTF, SLEVE, SNAP_NEAREST, or SNAP_BELOW. Not "
            << coordSysName << exit(FatalError);
}

/* For each column of fixed x value, find the point whose z value is closest to h.
   Update those points to have a z value of h.
 */
void snapNearestPointsToSurface(
        IOField<point>& newPoints,
        const Mountain& mountain)
{
	HashTable<int, scalar> minDistances;
	HashTable<int, scalar> closestZcoords;

	forAll(newPoints, ip)
	{
	    int x = roundUp(newPoints[ip].x(), 10); // FIXME: this hashtable method is dubious because it's not safe to compare doubles for equality
	    scalar z = newPoints[ip].z();
	    scalar h = mountain.heightAt(x);
	    scalar distance = abs(z - h);
	    
	    if (!minDistances.found(x) || distance < minDistances[x])
	    {
            minDistances.set(x, distance);
            closestZcoords.set(x, z);
	    }
	}

	forAll(newPoints, ip)
	{
	    int x = roundUp(newPoints[ip].x(), 10);
	    scalar z = newPoints[ip].z();
	    if (closestZcoords[x] == z)
	    {
		    newPoints[ip].z() = mountain.heightAt(x);
	    }
	}
}

void snapPointsBelowSurface(
        IOField<point>& newPoints,
        const Mountain& mountain,
        const scalar dz)
{
	forAll(newPoints, pointIdx)
	{
        scalar x = newPoints[pointIdx].x();
        scalar z = newPoints[pointIdx].z();
        scalar h = mountain.heightAt(x);
        if (z < h && z > h-dz)
        {
            newPoints[pointIdx].z() = h;
        }
        else if (z < h && z > h-2*dz)
        {
            newPoints[pointIdx].z() = h-0.01*dz;
        }
    }
}

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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
        mountainName == "AgnessiWitch" ? &AgnessiWitch :
        mountainName == "BottaKlein"   ? &BottaKlein : NULL;
        
    scalar (*fineMountain)(scalar, scalar) = 
        mountainName == "ScharExp" ? &ScharCos :
        mountainName == "ScharCos" ? &ScharCos :
        mountainName == "AgnessiWitch" ? &flatMountain :
        mountainName == "BottaKlein" ? &flatMountain : NULL;

    if 
    (
        mountainName != "ScharExp" && mountainName != "ScharCos"
     && mountainName != "AgnessiWitch" && mountainName != "BottaKlein"
    )
    {
        FatalErrorIn("add2dMountain")
      <<"mountainType should be one of ScharExp, ScharCos, BottaKlein or AgnessiWitch"
         << " not " << mountainName << exit(FatalError);
    }

    // Get which coord system to use
    enum coordSysType{BTF, HTF, SLEVE, SNAP_NEAREST, SNAP_BELOW, NONE};
    const word coordSysName(initDict.lookup("coordSys"));
    const coordSysType coordSys = coordSysName == "BTF" ? BTF :
                                  coordSysName == "HTF" ? HTF :
                                  coordSysName == "SLEVE" ? SLEVE : 
                                  coordSysName == "SNAP_NEAREST" ? SNAP_NEAREST : 
                                  coordSysName == "SNAP_BELOW" ? SNAP_BELOW : NONE;
    if (coordSys == NONE) badCoordinateSystem(coordSysName);
    
    // Declare and read in constants
    const scalar zt(readScalar(initDict.lookup("zt")));
    const scalar a(readScalar(initDict.lookup("a")));
    const scalar hm(readScalar(initDict.lookup("hm")));
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

    Mountain mountain(zt, a, hm, lambda, smoothMountain, fineMountain);
    
    // Calculate new points
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

    case SNAP_NEAREST:
        snapNearestPointsToSurface(newPoints, mountain);
        break;

    case SNAP_BELOW:
        {
            const scalar dz(readScalar(initDict.lookup("SNAP_dz")));
            snapPointsBelowSurface(newPoints, mountain, dz);
        }
        break;
    
    default:
        badCoordinateSystem(coordSysName);
    }

    newPoints.write();
}
