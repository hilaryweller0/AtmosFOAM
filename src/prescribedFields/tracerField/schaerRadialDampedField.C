#include "schaerRadialDampedField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(schaerRadialDampedField, 0);
addToRunTimeSelectionTable(tracerField, schaerRadialDampedField, dict);

schaerRadialDampedField::schaerRadialDampedField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
profileType(dict.lookupOrDefault<string>("profileType", string("default"))),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(0.001))),
rhoAir(dict.lookupOrDefault<scalar>("rhoAir", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3))),
p1(dict.lookupOrDefault<point>("ground", point(0, 0, -1000))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 1, 3e3)))
{};

scalar schaerRadialDampedField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    vector e = p - p1;
    scalar r = mag(vector(d.x()/A.x(), 0, d.z()/A.z()));

    if (r <= 1)
    {
        //If using specific ratio rather than relative ratio.
        if (profileType == "default")
        {
            return rho0*sqr(Foam::cos(M_PI*r/2)) + rhoAir;
        }
        else if (profileType == "rvs")
        {
            scalar theta0 = 300;
            scalar g = 9.81;
            scalar c_p = 1004;
            scalar muv = 0.018;
            scalar mud = 0.029;
    
            scalar T = 300*Foam::exp(-g*e.z()/(c_p*theta0));
            scalar P = 100000*Foam::exp(-g*e.z()/(c_p*theta0));
            scalar es = 611.2*Foam::exp( 17.67*(T-273.15)/(T-29.65) );
            
            return rho0 * muv/mud * es/(P-es) + rhoAir;
            //return rho0/rhoAir*sqr(Foam::cos(M_PI*r/2));
        }
        else if ((profileType == "rv") or (profileType == "rl"))
        {
            scalar theta0 = 300;
            scalar g = 9.81;
            scalar c_p = 1004;
            scalar muv = 0.018;
            scalar mud = 0.029;
    
            scalar T = 300*Foam::exp(-g*e.z()/(c_p*theta0));
            scalar P = 100000*Foam::exp(-g*e.z()/(c_p*theta0));
            scalar es = 611.2*Foam::exp( 17.67*(T-273.15)/(T-29.65) );
            
            if (profileType == "rv")
            {
                return min(muv/mud * es/(P-es), rho0/rhoAir*sqr(Foam::cos(M_PI*r/2)));
            }
            else if (profileType == "rl")
            {
                return rho0/rhoAir*sqr(Foam::cos(M_PI*r/2)) - min(muv/mud * es/(P-es), rho0/rhoAir*sqr(Foam::cos(M_PI*r/2))) ;
            }
        }
        else if (profileType == "q")
        {
            return rho0*sqr(Foam::cos(M_PI*r/2))/(rhoAir + rho0*sqr(Foam::cos(M_PI*r/2)));
        }
        else if (profileType == "block")
        {
            return rhoAir;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if ((profileType == "rv") or (profileType == "rl"))
        {
            return 0;
        }
        else
        {
            return rhoAir;
        }
    }
}
