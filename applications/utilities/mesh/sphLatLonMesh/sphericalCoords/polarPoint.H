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

C++
    listOfSourceFiles
        polarPoint.C
    endListOfSourceFiles

AUTHOR
    Hilary Spencer.
*/
// ------------------------------------------------------------------------- //

#ifndef polarPoint_H
#define polarPoint_H

#include "point.H"


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
public:

    // Static data members

        static scalar degToRad;
        static scalar radToDeg;
        
    // Constructors

        //- Construct null
        inline polarPoint() : VectorSpace<polarPoint, scalar, 3>(){}

        //- Construct from components
        inline polarPoint(const scalar lon_,const scalar lat_,const scalar r_);
        
        //- Construct from Istream
        inline polarPoint(Istream& is):VectorSpace<polarPoint, scalar, 3>(is){}

        //- Construct from vector
        inline polarPoint(const vector& v);


    // Member Functions


        // Access to components

        inline const scalar& lon() const { return component(0); }
        inline scalar& lon() { return component(0); }

        inline const scalar& lat() const { return component(1); }
        inline scalar& lat() { return component(1); }

        inline const scalar& r() const { return component(2); }
        inline scalar& r() { return component(2); }


        //- return the cartesian point
        inline point cartesian() const;

        //- return the cartesian point, given the global x, y and z directions
        inline point cartesian
        (
            const vector& xdir, const vector& ydir, const vector& zdir
        ) const;

        //- return the cartesian co-ordinates
        inline scalar x() const { return r()*cos(lon())*cos(lat()); }
        inline scalar y() const { return r()*sin(lon())*cos(lat()); }
        inline scalar z() const { return r()*sin(lat()); }

};

// non-member functions

//- convert cartesian point to polarPoint
inline polarPoint convertToPolar(const point& pc, const scalar maxLon = 360.);

//- convert polarPoint to cartesian point
inline point convertToCart(const polarPoint& pp);

//- convert cartesian point to polarPoint
inline polarPoint convertToPolar
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar maxLon = 360.
);

//- convert polarPoint to cartesian point
inline point convertToCart
(
    const scalar lon,
    const scalar lat,
    const scalar r
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polarPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
