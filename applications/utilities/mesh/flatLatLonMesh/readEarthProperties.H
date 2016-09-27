Info << "\nReading earthProperties" << endl;

IOdictionary earthProperties
(
    IOobject
    (
        "earthProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

//dimensionedScalar earthRadius(earthProperties.lookup("earthRadius"));
////dimensionedScalar atmosHeight(earthProperties.lookup("atmosHeight"));

//const label nLat(readLabel(earthProperties.lookup("nLat")));
//const label nLon(readLabel(earthProperties.lookup("nLon")));

//const scalar maxLonRatio = earthProperties.found("maxLonRatio") ?
//                      readScalar(earthProperties.lookup("maxLonRatio")) : 0;

//scalarList urefLats = earthProperties.found("urefLats") ?
//              scalarList(earthProperties.lookup("urefLats")) : scalarList(0);

//const label nLayers = earthProperties.found("nLayers") ?
//        readLabel(earthProperties.lookup("nLayers")) : 1;

//const vector meshRotation = earthProperties.found("meshRotation") ?
//  vector(earthProperties.lookup("meshRotation"))*constant::mathematical::pi/180. :
//  vector::zero;

//Info << "Creating spherical polar mesh with " << nLat << " latitudes, " << nLon
//<< " longitudes, maxLonRatio = " << maxLonRatio << " urefLats = " << urefLats
//<< " nLayers = " << nLayers << " and meshRotation = " << meshRotation << endl;
