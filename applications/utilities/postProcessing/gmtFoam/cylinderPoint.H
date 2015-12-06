// The FOAM Project // File: cylinderPoint.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   cylinderPoint
   \\  /           | Family: 
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    cylinderPoint

DESCRIPTION
    Stores a point in cylindrical polar co-ordinates as theta, r, z

*/
// ------------------------------------------------------------------------- //

#ifndef cylinderPoint_H
#define cylinderPoint_H

#include "point.H"
#include "mathematicalConstants.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class cylinderPoint Declaration
\*---------------------------------------------------------------------------*/

class cylinderPoint
:
    public VectorSpace<cylinderPoint, scalar, 3>
{
public:

    // Static data members

        static scalar const degToRad = M_PI/180.;
        static scalar const radToDeg = 180./M_PI;
        
    // Constructors

        //- Construct null
        inline cylinderPoint():VectorSpace<cylinderPoint, scalar, 3>(){}

        //- Construct given cylinderPoint
        inline cylinderPoint(const cylinderPoint& pp)
            :VectorSpace<cylinderPoint, scalar, 3>(pp){}

        //- Construct given three components
        inline cylinderPoint(const scalar theta_,const scalar r_,const scalar z_)
            :VectorSpace<cylinderPoint, scalar, 3>()
        {
            theta() = theta_;
            r() = r_;
            z() = z_;
        }

        //- Construct from Istream
        inline cylinderPoint(Istream& is)
            :VectorSpace<cylinderPoint, scalar, 3>(is){}

    // Member Functions

        // Access to components

        inline const scalar& theta() const { return component(0); }
        inline scalar& theta() { return component(0); }

        inline const scalar& r() const { return component(1); }
        inline scalar& r() { return component(1); }

        inline const scalar& z() const { return component(2); }
        inline scalar& z() { return component(2); }


        //- return the cartesian point
        inline point cartesian() const
        {
            return point
            (
                r()*cos(theta()),
                r()*sin(theta()),
                z()
            );
        }

        //- return the cartesian point, given the global x, y and z directions
        inline point cartesian
        (
            const vector& xdir, const vector& ydir, const vector& zdir
        ) const
        {
            point p
            (
                r()*cos(theta()),
                r()*sin(theta()),
                z()
            );
            
            tensor C(xdir, ydir, zdir);
            tensor invC = mag(det(C)) > VSMALL ? inv(C) : tensor::zero;

            return point(invC & p);
        }

        //- return the cartesian co-ordinates
        inline scalar x() const { return r()*cos(theta()); }
        inline scalar y() const { return r()*sin(theta()); }

};

// non-member functions

//- convert cartesian point to cylinderPoint
inline cylinderPoint convertToCylinder
(
    const point& pc, const scalar maxTheta = 360.
)
{
    using namespace constant::mathematical;

    const scalar maxThetaR = maxTheta*cylinderPoint::degToRad;
    cylinderPoint pp;
    
    pp.z() = pc.z();
    pp.r() = mag(vector(pc.x(),pc.y(),0));
    pp.theta() = atan2(pc.y(), pc.x());

    while(pp.theta() > maxThetaR) {pp.theta() -= 2*pi;}
    while(pp.theta() < maxThetaR - 2*pi) {pp.theta() += 2*pi;}

    return pp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
