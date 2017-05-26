#include "temperatureGradientField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(temperatureGradientField, 0);
addToRunTimeSelectionTable(tracerField, temperatureGradientField, dict);

temperatureGradientField::temperatureGradientField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
tracerField(velocityField),
rho0(dict.lookupOrDefault<scalar>("maxMagnitude", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(-50e3, 0, 9e3))),
A(dict.lookupOrDefault<vector>("halfWidth", vector(25e3, 1, 3e3)))
{};

scalar temperatureGradientField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    scalar theta0 = 300;
    scalar g = 9.81;
    scalar c_p = 1004;
    
    return theta0*Foam::exp(-g*d.z()/(c_p*theta0));

}
