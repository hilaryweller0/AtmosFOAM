// The FOAM Project // File: polarVector.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   polarVector
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    polarVector

DESCRIPTION
    stores a vector in global spherical polar co-ordinates 
    as lonDot, latDot, rDot
    
C++
    listOfSourceFiles
        polarVector.C
    endListOfSourceFiles

AUTHOR
    Hilary Spencer.
*/
// ------------------------------------------------------------------------- //

#ifndef polarVector_H
#define polarVector_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class forward declaration
class lCartPolarVector;
class polarPoint;

/*---------------------------------------------------------------------------*\
                           Class polarVector Declaration
\*---------------------------------------------------------------------------*/

class polarVector
:
    public VectorSpace<polarVector, scalar, 3>
{

public:

    // Constructors

        //- Construct null
        polarVector() : VectorSpace<polarVector, scalar, 3>(){}

        //- Construct from components
        polarVector
        (
            const scalar lonDot,
            const scalar latDot, 
            const scalar rDot
        );
        
        //- Construct from Istream
        polarVector(Istream& is) : VectorSpace<polarVector, scalar, 3>(is){}


    // Member Functions

        // Access to components
        
        const scalar& lonDot() const { return component(0); }
        scalar& lonDot() { return component(0); }
        
        const scalar& latDot() const { return component(1); }
        scalar& latDot() { return component(1); }
        
        const scalar& rDot() const { return component(2); }
        scalar& rDot() { return component(2); }
};

// * * * * * * * * non-member conversion functions * * * * * * * * * * * * * //

// convert to a polarVector

    //- convert cartesian vector and cartesian point to polarVector
    //polarVector convertToPolar(const vector& vc, const point& pc);

    //- convert local carteisan polar vector and cartesian point to polarVector
    polarVector convertToPolar(const lCartPolarVector& lcv, const point& pc);

    //- convert cartesian vector and polarPoint to polarVector
    //polarVector convertToPolar(const vector& vc, const polarPoint& pp);

    //- convert local carteisan polar vector and polarPoint to polarVector
    polarVector convertToPolar
    (
        const lCartPolarVector& lcv, 
        const polarPoint& pp
    );


// convert to a cartesian vector

    //- convert polarVector at a carteisan point to a cartesian vector
    vector convertToCart(const polarVector& pv, const point& pc);

    //- convert polarVector at a polarPoint to a cartesian vector
    vector convertToCart(const polarVector& pv, const polarPoint& pp);

    //- convert local cartesian vector at a carteisan point to cartesian vector
    vector convertToCart(const lCartPolarVector& lcv, const point& pc);

    //- convert local cartesian vector at a polarPoint to a cartesian vector
    vector convertToCart(const lCartPolarVector& lcv, const polarPoint& pp);


// convert to a local cartesian vector

    //- convert polarVector at a carteisan point to a local cartesian vector
    lCartPolarVector convertToLCart(const polarVector& pv, const point& pc);

    //- convert polarVector at a polarPoint to a local cartesian vector
    lCartPolarVector convertToLCart
    (
        const polarVector& pv, 
        const polarPoint& pp
    );

    //- convert vector at a carteisan point to a local cartesian vector
    lCartPolarVector convertToLCart(const vector& vc, const point& pc);

    //- convert vector at a polarPoint to a local cartesian vector
    lCartPolarVector convertToLCart
    (
        const vector& vc, 
        const polarPoint& pp
    );


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
