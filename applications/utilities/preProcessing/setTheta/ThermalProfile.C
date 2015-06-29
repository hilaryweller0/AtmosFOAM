#include "ThermalProfile.H"
#include "fvCFD.H"

ThermalProfile::ThermalProfile(
            const IOdictionary& environmentalProperties,
            const dimensionedVector g,
            const dimensionedScalar T0) :
    g(g),
    T0(T0),
    nLayers(environmentalProperties.lookup("BruntVaisallaFreq")),
    zLayers(environmentalProperties.lookup("zN"))
{
    if (nLayers.size()+1 != zLayers.size())
    {
        FatalErrorIn("ThermalProfile")
            << " size of BruntVaisallaFreq in environmentalProperties should be"
            << " one smaller than the size of zN"
            << exit(FatalError);
    }
}

scalar ThermalProfile::thetaAt(const point& p) const
{
    const scalar z = p.z();

    scalar theta0 = T0.value();
    for(label il = 0; il < nLayers.size(); il++)
    {
        if (z >= zLayers[il] && z < zLayers[il+1])
        {
            return theta0*Foam::exp
            (
                sqr(nLayers[il])/mag(g.value())*(z-zLayers[il])
            );
        }
        theta0 *= Foam::exp
            (sqr(nLayers[il])/mag(g.value())*(zLayers[il+1]-zLayers[il]));
    }
    FatalErrorIn("ThermalProfile") << "height " << z << " not found in levels "
        << zLayers << exit(FatalError);
    return -1;
}

