// The FOAM Project // File: polarPoint.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   polarPoint
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    polarPoint

DESCRIPTION
    stores a point in global spherical polar co-ordinates as lon, lat, r

*/
// ------------------------------------------------------------------------- //

#ifndef polarPoint_H
#define polarPoint_H

#include "point.H"
#include "mathematicalConstants.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class polarPoint Declaration
\*---------------------------------------------------------------------------*/

class polarPoint
:
    public VectorSpace<polarPoint, scalar, 3>
{
    // Private data
        // Earth radius
        scalar a_;
public:

    // Static data members

        static scalar constexpr degToRad = M_PI/180.;
        static scalar constexpr radToDeg = 180./M_PI;
        
    // Constructors

        //- Construct null
        inline polarPoint():VectorSpace<polarPoint, scalar, 3>(),a_(0){}

        //- Construct given polarPoint
        inline polarPoint(const polarPoint& pp)
            :VectorSpace<polarPoint, scalar, 3>(pp),a_(pp.a_){}

        //- Construct given three (or 4) components
        inline polarPoint(const scalar lon_,const scalar lat_, const scalar r_,
                          const scalar a__ = 0)
            :VectorSpace<polarPoint, scalar, 3>(),
             a_(a__)
        {
            lon() = lon_;
            lat() = lat_;
            r() = r_;
        }

        //- Construct from Istream
        inline polarPoint(Istream& is)
            :VectorSpace<polarPoint, scalar, 3>(is), a_(readScalar(is)){}

    // Member Functions

        // Access to components

        inline const scalar& lon() const { return component(0); }
        inline scalar& lon() { return component(0); }

        inline const scalar& lat() const { return component(1); }
        inline scalar& lat() { return component(1); }

        inline const scalar& r() const { return component(2); }
        inline scalar& r() { return component(2); }

        inline const scalar& a() const { return a_; }
        inline scalar& a() { return a_; }

        inline scalar Lon() const { return lon()*radToDeg; }
        inline scalar Lat() const { return lat()*radToDeg; }
        inline scalar Z() const { return r() - a(); }

        //- return the cartesian point
        inline point cartesian() const
        {
            const scalar coslat = cos(lat());

            return point
            (
                r()*cos(lon())*coslat,
                r()*sin(lon())*coslat,
                r()*sin(lat())
            );
        }

        //- return the cartesian point, given the global x, y and z directions
        inline point cartesian
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

        //- return the cartesian co-ordinates
        inline scalar x() const { return r()*cos(lon())*cos(lat()); }
        inline scalar y() const { return r()*sin(lon())*cos(lat()); }
        inline scalar z() const { return r()*sin(lat()); }

};

// non-member functions

//- convert cartesian point to polarPoint
inline polarPoint convertToPolar
(
    const point& pc,
    const scalar maxLon = 360.,
    const scalar a__ = 0
)
{
    using namespace constant::mathematical;

    const scalar maxLonR = maxLon*polarPoint::degToRad;
    polarPoint pp;
    pp.r() = mag(pc);

    scalar sinlat = max(min(pc.z() / pp.r(), 1.), -1.);
    pp.lat() = asin(sinlat);

    pp.lon() = atan2(pc.y(), pc.x());

    while(pp.lon() > maxLonR) {pp.lon() -= 2*pi;}
    while(pp.lon() < maxLonR - 2*pi) {pp.lon() += 2*pi;}
    
    pp.a() = a__;

    return pp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
