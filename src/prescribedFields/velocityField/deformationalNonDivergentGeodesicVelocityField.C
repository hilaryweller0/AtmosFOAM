#include "deformationalNonDivergentGeodesicVelocityField.H"
#include "addToRunTimeSelectionTable.H"
#include "polarPoint.H"
#include "deformationalNonDivergentGeodesicVelocityFieldData.H"

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
    radius(readScalar(dict.lookup("radius"))),
    deformationScale(readScalar(dict.lookup("deformationScale"))),
    endTime(readScalar(dict.lookup("endTime")))
{};

vector deformationalNonDivergentGeodesicVelocityField::streamfunctionAt
(
    const label ip,
    const sphericalMeshData& spherical,
    const Time& t
) const
{
    // Get reference to the data
    const deformationalNonDivergentGeodesicVelocityFieldData&
         data = deformationalNonDivergentGeodesicVelocityFieldData::New
    (
        spherical.mesh(), radius, deformationScale, endTime
    );

    const scalar lon = spherical.pointsLatLonz()[ip][0];

    // section 2.3 doi:10.5194/gmd-5-887-2012
    const scalar lonPrime = lon - 2 * M_PI * t.value() / endTime;
    
    return data.phat_sqrrByT()[ip]*
    (
        data.Dcos2Lat()[ip]*sqr(Foam::sin(lonPrime))
            *Foam::cos(M_PI*t.value()/endTime)
      - data.TwoPiSinLat()[ip]
      + 2 * M_PI
    );
}
