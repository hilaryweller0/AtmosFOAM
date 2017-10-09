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
T0(dict.lookupOrDefault<scalar>("centralValue", scalar(1))),
p0(dict.lookupOrDefault<point>("centre", point(0, 0, 0))),
gradientType(dict.lookupOrDefault<string>("gradientType", string("exponential"))),
g(dict.lookupOrDefault<scalar>("gravitationalConstant", scalar(9.81))),
cp(dict.lookupOrDefault<scalar>("specificHeatConstPressure", scalar(1004))),
theta0(dict.lookupOrDefault<scalar>("decayConstant", scalar(300)))
{};

scalar temperatureGradientField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    vector d = p - p0;
    
    if (gradientType == "exponential")
    {
        return T0*Foam::exp(-g*d.z()/(cp*theta0));
    }
    else if (gradientType == "linear")
    {
        return T0 - g*d.z()/(cp*theta0);
    }
    else
    {
        return T0;
    }
}
