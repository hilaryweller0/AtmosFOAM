#include "deformationalNonDivergentGeodesicVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "polarPoint.H"
#include "deformationalNonDivergentGeodesicVelocityFieldData.H"

namespace Foam
{
defineTypeNameAndDebug(deformationalNonDivergentGeodesicVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    deformationalNonDivergentGeodesicVelocityField,
    dict
);

deformationalNonDivergentGeodesicVelocityField::
deformationalNonDivergentGeodesicVelocityField(const dictionary& dict)
:
    geodesicNonDivergentVelocityField(dict),
    deformationScale(readScalar(dict.lookup("deformationScale")))
{};

vector deformationalNonDivergentGeodesicVelocityField::streamfunctionAt
(
    const label ip,
    const sphericalMeshData& spherical,
    scalar time
) const
{
    // Get reference to the data
    const deformationalNonDivergentGeodesicVelocityFieldData&
         data = deformationalNonDivergentGeodesicVelocityFieldData::New
    (
        spherical.mesh(), earthRadius_.value(), deformationScale, endTime_.value()
    );

    const scalar lon = spherical.pointsLatLonz()[ip][0];

    // section 2.3 doi:10.5194/gmd-5-887-2012
    const scalar lonPrime = lon - 2 * M_PI * time / endTime_.value();
    
    return data.phat_sqrrByT()[ip]*
    (
        data.Dcos2Lat()[ip]*sqr(Foam::sin(lonPrime))
            *Foam::cos(M_PI*time/endTime_.value())
      - data.TwoPiSinLat()[ip]
      + 4 * M_PI
    );
}
}
