// Check to see if the projection is G or g (Orthographic) and if so set the
// limits lat and lon are written out for

const bool ortho = 
(
    projection[0] == 'G' || projection[0] == 'g'
// || projection[0] == 'S' || projection[0] == 's'
// || projection[0] == 'Q' || projection[0] == 'q'
// || projection[0] == 'W' || projection[0] == 'w'
) ? true : false;
scalar lon0 = 0;
scalar lat0 = 0;

if (ortho)
{
    string lonLatStr = projection;
    lonLatStr.replace(projection[0], "");
    lonLatStr.replaceAll("/", " ");
    IStringStream lonLat(lonLatStr);
    lonLat >> lon0 >> lat0;
    lon0 /= radToDeg;
    lat0 /= radToDeg;
}

//const bool latLon = projection2D == "latLon";

// If projection includes latitute or longitude, read in Omega to define
// latitude and longitude
// The axis of rotation
vector axis = vector(0,0,1);
vector xdir = vector::zero;
vector ydir = vector::zero;
if (projection2D == "latLon" || projection2D == "latz" || projection2D == "lonz")
{
    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dimensionedVector Omega
    (
        earthProperties.lookupOrDefault<dimensionedVector>
        (
            "Omega",
            dimensionedVector("", dimless, vector(0,0,1))
        )
    );

    axis = Omega.value();
    if (mag(axis) > SMALL)
    {
        axis /= mag(axis);
    }
    else axis = vector(0,0,1);
    
    if (mag(axis - vector(1,0,0)) > SMALL)
    {
        xdir = vector(1,0,0);
    }
    else xdir = vector(0,1,0);
    
    xdir -= (xdir&axis)*axis;
    xdir /= mag(xdir);
    
    ydir = axis ^ xdir;
    
    Info << "xdir = " << xdir << " ydir = " << ydir << " axis = " << axis
         << endl;
}

