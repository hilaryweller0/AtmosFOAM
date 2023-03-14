#include "geodesicSlottedCylinderTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(geodesicSlottedCylinderTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicSlottedCylinderTracerField, dict);

geodesicSlottedCylinderTracerField::geodesicSlottedCylinderTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    Rsphere(readScalar(dict.lookup("Rsphere"))),
    Rcylinder(readScalar(dict.lookup("Rcylinder"))),
    hBackground(dict.lookupOrDefault<scalar>("hBackground", scalar(0))),
    hmax(readScalar(dict.lookup("hmax")))
{};

scalar geodesicSlottedCylinderTracerField::tracerAt
(
    const point& p,
    const Time& t
) const
{
    scalar lon1 = 5.0*M_PI/6.0;
    scalar lat1 = 0;

    scalar lon2 = -5.*M_PI/6.0;
    scalar lat2 = 0;
    
    scalar plat = Foam::asin(p.z()/mag(p));
    scalar plon = Foam::atan2(p.y(), p.x());

    const point centre1
    (
        Rsphere * Foam::cos(lat1) * Foam::cos(lon1),
        Rsphere * Foam::cos(lat1) * Foam::sin(lon1),
        Rsphere * Foam::sin(lat1)
    );

    const point centre2
    (
        Rsphere * Foam::cos(lat2) * Foam::cos(lon2),
        Rsphere * Foam::cos(lat2) * Foam::sin(lon2),
        Rsphere * Foam::sin(lat2)
    );
    
    scalar tracer = hBackground;
    scalar r1 = mag(p - centre1); 
    scalar r2 = mag(p - centre2);
    
    if
    (
        (r1 <= Rcylinder && mag(lon1 - plon) >= Rcylinder/(6*Rsphere))
     || (r2 <= Rcylinder && mag(lon2 - plon) >= Rcylinder/(6*Rsphere))
     || (r1 <= Rcylinder && mag(lon1 - plon) < Rcylinder/(6*Rsphere) && plat - lat1 < -5./12*Rcylinder/Rsphere)
     || (r2 <= Rcylinder && mag(lon2 - plon) < Rcylinder/(6*Rsphere) && plat - lat2 > 5./12*Rcylinder/Rsphere)
    )
    {
        tracer = hmax;
    }
    
    return tracer;
}
