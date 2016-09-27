// The FOAM Project // File: polarVector.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   polarVector
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    stores a vector in global spherical polar co-ordinates 
    as lonDot, latDot, rDot

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/

#include "polarVector.H"
#include "polarPoint.H"
#include "lCartPolarVector.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polarVector::polarVector
(
    const scalar lonDot,
    const scalar latDot, 
    const scalar rDot
)
{
    component(0) = lonDot;
    component(1) = latDot;
    component(2) = rDot;
}        


// * * * * * * * * * * non-Member conversion Functions  * * * * * * * * * //

// * * * * * * * * * * convert to a polarVector  * * * * * * * * * //

// convert cartesian vector and cartesian point to polarVector
//Foam::polarVector Foam::convertToPolar(const vector& vc, const point& pc)
//{}

// convert cartesian vector and polarPoint to polarVector
//Foam::polarVector Foam::convertToPolar(const vector& vc, const polarPoint& pp)
//{}

// convert local carteisan polar vector and cartesian point to polarVector
Foam::polarVector Foam::convertToPolar
(
    const lCartPolarVector& lcv, 
    const point& pc
)
{
    return convertToPolar(lcv, convertToPolar(pc));
}

// convert local carteisan polar vector and polarPoint to polarVector
Foam::polarVector Foam::convertToPolar
(
    const lCartPolarVector& lcv, 
    const polarPoint& pp
)
{
    const scalar londot = lcv.us()/(pp.r()*cos(pp.lat()));
    const scalar latdot = lcv.vs()/pp.r();
    
    return polarVector(londot, latdot, lcv.ws());
}

// * * * * * * * * * * convert to a cartesian vector  * * * * * * * * * //

// convert polarVector at a carteisan point to a cartesian vector
Foam::vector Foam::convertToCart(const polarVector& pv, const point& pc)
{
    return convertToCart(pv, convertToPolar(pc));
}


// convert polarVector at a polarPoint to a cartesian vector
Foam::vector Foam::convertToCart(const polarVector& pv, const polarPoint& pp)
{
    const scalar coslat = cos(pp.lat());
    const scalar coslon = cos(pp.lon());
    const scalar sinlat = sin(pp.lat());
    const scalar sinlon = sin(pp.lon());
    
    const scalar uc = pv.rDot()*coslon*coslat 
                     - pp.r()*pv.lonDot()*sinlon*coslat
                     - pp.r()*pv.latDot()*coslon*sinlat;

    const scalar vc = pv.rDot()*sinlon*coslat 
                     + pp.r()*pv.lonDot()*coslon*coslat
                     - pp.r()*pv.latDot()*sinlon*sinlat;

    const scalar wc = pv.rDot()*sinlat + pp.r()*pv.latDot()*coslat;
    
    return vector(uc, vc, wc);
}

// convert local cartesian vector at a carteisan point to cartesian vector
Foam::vector Foam::convertToCart(const lCartPolarVector& lcv, const point& pc)
{
    return convertToCart(lcv, convertToPolar(pc));
}

// convert local cartesian vector at a polarPoint to a cartesian vector
Foam::vector Foam::convertToCart
(
    const lCartPolarVector& lcv, 
    const polarPoint& pp
)
{
    const scalar coslat = cos(pp.lat());
    const scalar coslon = cos(pp.lon());
    const scalar sinlat = sin(pp.lat());
    const scalar sinlon = sin(pp.lon());
    
    const scalar londot = lcv.us()/(pp.r()*coslat);
    const scalar latdot = lcv.vs()/pp.r();
    const scalar rdot = lcv.ws();
    
    const scalar uc = rdot*coslon*coslat - pp.r()*londot*sinlon*coslat
                                         - pp.r()*latdot*coslon*sinlat;
    const scalar vc = rdot*sinlon*coslat + pp.r()*londot*coslon*coslat
                                         - pp.r()*latdot*sinlon*sinlat;
    const scalar wc = rdot*sinlat + pp.r()*latdot*coslat;
    
    return vector(uc, vc, wc);
}


// * * * * * * * * * * convert to a local cartesian vector  * * * * * * * * //

// convert polarVector at a carteisan point to a local cartesian vector
Foam::lCartPolarVector Foam::convertToLCart
(
    const polarVector& pv, 
    const point& pc
)
{
    return convertToLCart(pv, convertToPolar(pc));
}


// convert polarVector at a polarPoint to a local cartesian vector
Foam::lCartPolarVector Foam::convertToLCart
(
    const polarVector& pv, 
    const polarPoint& pp
)
{
    return lCartPolarVector
    (
        pv.lonDot() * pp.r() * cos(pp.lat()), 
        pv.latDot() * pp.r(), 
        pv.rDot()
    );
}

// convert vector at a carteisan point to a local cartesian vector
Foam::lCartPolarVector Foam::convertToLCart
(
    const vector& vc,
    const point& pc
)
{
    return convertToLCart(vc, convertToPolar(pc));
}


// convert vector at a polarPoint to a local cartesian vector
Foam::lCartPolarVector Foam::convertToLCart
(
    const vector& vc, 
    const polarPoint& pp
)
{
    return lCartPolarVector
    (
        -vc[0]*sin(pp.lon()) + vc[1]*cos(pp.lon()),
        -vc[0]*cos(pp.lon())*sin(pp.lat()) - vc[1]*sin(pp.lon())*sin(pp.lat())
                                           + vc[2]*cos(pp.r()),
        -vc[0]*cos(pp.lon())*cos(pp.lat()) - vc[1]*sin(pp.lon())*cos(pp.lat())
                                           + vc[2]*sin(pp.r())
    );
}

// ************************************************************************* //
