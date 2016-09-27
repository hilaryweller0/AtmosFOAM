// The FOAM Project // File: polarPointI.H
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   polarPoint
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/

#include "mathematicalConstants.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline Foam::polarPoint::polarPoint
(
    const scalar lon_,
    const scalar lat_,
    const scalar r_
)
{
    lon() = lon_;
    lat() = lat_;
    r() = r_;
}

// Construct from vector
inline Foam::polarPoint::polarPoint(const vector& v)
{
    lon() = v.component(0);
    lat() = v.component(1);
    r() = v.component(2);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// return the cartesian point
inline Foam::point Foam::polarPoint::cartesian() const
{
    const scalar coslat = cos(lat());

    return point
    (
        r()*cos(lon())*coslat,
        r()*sin(lon())*coslat,
        r()*sin(lat())
    );
}

// return the cartesian point given the global x, y and z directions
inline Foam::point Foam::polarPoint::cartesian
(
    const vector& xdir, const vector& ydir, const vector& zdir
) const
{
    const scalar coslat = cos(lat());

    point p
    (
        r()*cos(lon())*coslat,
        r()*sin(lon())*coslat,
        r()*sin(lat())
    );
    
    tensor C(xdir, ydir, zdir);
    tensor invC = mag(det(C)) > VSMALL ? inv(C) : tensor::zero;

    return point(invC & p);
}

// * * * * * * * * * * * * * * * non-Member Functions  * * * * * * * * * * * //

// convert cartesian point to polarPoint
inline Foam::polarPoint Foam::convertToPolar
(
    const point& pc, const scalar maxLon
)
{
    using namespace constant::mathematical;

    const scalar maxLonR = maxLon*polarPoint::degToRad;

    polarPoint pp;

    pp.r() = mag(pc);

    scalar sinlat = max(min(pc.z() / pp.r(), 1.), -1.);
    pp.lat() = asin(sinlat);

    pp.lon() = atan2(pc.y(), pc.x());

//    scalar yOverRcoslat = pc.y() / ( pp.r() * cos(pp.lat()));

//    if(yOverRcoslat >= 1)       pp.lon() = 0.5*pi;
//    else if(yOverRcoslat <= -1) pp.lon() = -0.5*pi;
//    else if(pc.x() > 0)         pp.lon() = asin(yOverRcoslat);
//    else                        pp.lon() = pi - asin(yOverRcoslat);

    while(pp.lon() > maxLonR) {pp.lon() -= 2*pi;}
    while(pp.lon() < maxLonR - 2*pi) {pp.lon() += 2*pi;}

    return pp;
}

// convert polarPoint to cartesian point
inline Foam::point Foam::convertToCart(const polarPoint& pp)
{
    const scalar coslat = cos(pp.lat());

    return point
    (
        pp.r()*cos(pp.lon())*coslat,
        pp.r()*sin(pp.lon())*coslat,
        pp.r()*sin(pp.lat())
    );
}


// convert cartesian point to polarPoint
inline Foam::polarPoint Foam::convertToPolar
(
    const scalar x, const scalar y, const scalar z, const scalar maxLon
)
{
    using namespace constant::mathematical;

    const scalar maxLonR = maxLon*polarPoint::degToRad;

    polarPoint pp;

    pp.r() = sqrt(sqr(x) + sqr(y) + sqr(z));

    scalar sinlat = max(min(z / pp.r(), 1.), -1.);
    pp.lat() = asin(sinlat);

    pp.lon() = atan2(y,x);

    while(pp.lon() > maxLonR) {pp.lon() -= 2*pi;}
    while(pp.lon() < maxLonR - 2*pi) {pp.lon() += 2*pi;}

    return pp;
}

// convert polarPoint to cartesian point
inline Foam::point Foam::convertToCart
(
    const scalar lon, const scalar lat, const scalar r
)
{
    const scalar coslat = cos(lat);

    return point
    (
        r*cos(lon)*coslat,
        r*sin(lon)*coslat,
        r*sin(lat)
    );
}



// ************************************************************************* //
