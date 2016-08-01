/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef VectorSpaceFunctions_H
#define VectorSpaceFunctions_H

#include "error.H"
#include "products.H"
#include "vector.H"
#include "ops.H"
#include "mathematicalConstants.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// The unit vector of a vector
template<class Form, class Cmpt, direction nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt> unitVector
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt m = mag(vs);
    if (m < VSMALL) return vs;
    
    return vs/m;
}

//- Approximate distance between 2 points on a sphere
inline scalar angleBetween(const point& a, const point& b)
{
    const scalar ra = mag(a);
    const scalar rb = mag(b);

    // Angle does not exixt if one point has zero magnitude
    if (ra < VSMALL || rb < VSMALL)
    {
        FatalErrorIn("angleBetween") << " arc length between points " << a << " and "
            << b << " does not exist because one is zero" << exit(FatalError);
    }
    
    // angle between the points
    const scalar theta = 2*asin(0.5*mag(a/ra - b/rb));
    
    return theta;
}

//- Approximate distance between 2 points on a sphere
inline scalar arcLength(const point& a, const point& b)
{
    const scalar ra = mag(a);
    const scalar rb = mag(b);
    const scalar rdiff = ra - rb;
    if (mag(rdiff) > ra*1e-9 || ra < VSMALL || rb < VSMALL)
    {
        FatalErrorIn("arcLength")
            <<" trying to find arc length between points with different radii"
            << setprecision(12)
            << "\nmag(a) = " << ra << " a = " << a
            << "\nmag(b) = " << rb << " b = " << b
            << "\ndiff = " << rdiff
            << exit(FatalError);
    }
    const point pa = a/ra;
    const point pb = b/rb;
    // angle between the points
    const scalar theta = 2*asin(0.5*mag(pa - pb));
    
    return ra*theta;
}

//- Distance between 2 points in Cartesian geometry
inline scalar dist(const point& a, const point& b)
{
    return mag(a - b);
}


//- Approximate distance between 2 points in spherical geometry
inline scalar sphDist(const point& a, const point& b)
{
    const scalar ra = mag(a);
    const scalar rb = mag(b);
    const scalar rdiff = ra - rb;
    const scalar rmid = 0.5*(ra + rb);
    const point pa = a/ra;
    const point pb = b/rb;
    // angle between the points
    const scalar theta = 2*asin(0.5*mag(pa - pb));
    
    return sqrt(sqr(rmid*theta) + sqr(rdiff));
}


inline scalar sphTriSolidAngle(const point& a, const point& b, const point& c)
{
    // pa, pb and pc are the points on the unit sphere
    const point pa = unitVector(a);
    const point pb = unitVector(b);
    const point pc = unitVector(c);
    // A, B and C are the great circle distances between points (or angles)
    // formula using sin rather than cos for less roundoff
    const scalar A = 2*asin(0.5*mag(pb - pc));
    const scalar B = 2*asin(0.5*mag(pa - pc));
    const scalar C = 2*asin(0.5*mag(pa - pb));
    // formular for the solid angle using the excess
    const scalar s = 0.5*(A + B + C);
    const scalar t = tan(0.5*s)*tan(0.5*(s-A))*tan(0.5*(s-B))*tan(0.5*(s-C));
    const scalar E = 4*atan(sqrt(mag(t)));
    return E;
}


inline scalar sphTriArea(const point& a, const point& b, const point& c)
{
    return magSqr(a)*sphTriSolidAngle(a, b, c);
}

inline scalar sphTriDistAngle(const point& a, const point& b, const point& c)
{
    // pa, pb and pc are the points on the unit sphere
    const point pa = unitVector(a);
    const point pb = unitVector(b);
    const point pc = unitVector(c);
    
    // vectors between points
    vector ab = pb - pa;
    vector ac = pc - pa;
    
    // Great circle distances between points
    const scalar AB = 2*asin(0.5*mag(ab));
    const scalar AC = 2*asin(0.5*mag(ac));
    
    // make vectors unit vectors
    ab = unitVector(ab);
    ac = unitVector(ac);

    // unit triangle area
    return 0.5*AB*AC*mag(ab ^ ac);
}

// exponential map to map a point xi on the surface of a sphere a distance
// |map| in direction map to another point on the sphere
inline vector expMapOnSphere(const vector& xi, const vector& map)
{
    scalar r = mag(xi);
    if (r < SMALL) return vector::zero;
    scalar theta = mag(map)/r;
    if (theta < SMALL) return xi;
    
    vector e = unitVector(xi ^ map);
    tensor transform = tensor::I*cos(theta)
                     + (1 - cos(theta))*(e*e)
                     + tensor(0,  -e[2],e[1],
                              e[2],0,  -e[0],
                             -e[1],e[0],0)*sin(theta);
    
    return transform & xi;
    
//    tensor transform(xi, map, (xi ^ map));
//    vector x = sqr(r)*(inv(transform)&vector(cos(theta), theta*sin(theta), 0));
//    vector check = (transform & x);
//    
//    Info << "\nMapping point " << xi << " with map " << map
//         << nl << " with transform " << transform << " r = " << r
//         << " theta = " << theta << nl << " to point "
//         << x << " with radius " << mag(x)
//         << nl << "det(transform) = " << det(transform)
//         << " check = " << check << " RHS = " 
//         << sqr(r)*vector(cos(theta), theta*sin(theta), 0)
//         << nl << "inv(transform) = " << inv(transform)
//         << nl << "inv(transform) * transform = "
//         << inv(transform) & transform << nl;
//    
//    return r*(inv(transform)&vector(cos(theta), theta*sin(theta), 0));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
