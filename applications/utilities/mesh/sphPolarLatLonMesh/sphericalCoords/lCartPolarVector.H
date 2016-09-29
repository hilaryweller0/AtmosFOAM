// The FOAM Project // File: lCartPolarVector.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   lCartPolarVector
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    lCartPolarVector

DESCRIPTION
    stores a vector in local cartesian co-ordinates on the sphere
    and us, vs, ws (being the velocity components in the lon, lat, r directions)

C++
    listOfSourceFiles
        lCartPolarVector.C
    endListOfSourceFiles

AUTHOR
    Hilary Spencer.
*/
// ------------------------------------------------------------------------- //

#ifndef lCartPolarVector_H
#define lCartPolarVector_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class forward declaration
class polarVector;
class polarPoint;

/*---------------------------------------------------------------------------*\
                           Class lCartPolarVector Declaration
\*---------------------------------------------------------------------------*/

class lCartPolarVector
:
    public VectorSpace<lCartPolarVector, scalar, 3>
{

public:

    // Constructors

        //- Construct null
        lCartPolarVector() : VectorSpace<lCartPolarVector, scalar, 3>(){}

        //- Construct from components
        lCartPolarVector
        (
            const scalar us,
            const scalar vs, 
            const scalar rs
        );
        
        //- Construct from Istream
        lCartPolarVector(Istream& is) 
            : VectorSpace<lCartPolarVector, scalar, 3>(is){}


    // Member Functions

        // Access to components
        
        const scalar& us() const { return component(0); }
        scalar& us() { return component(0); }
        
        const scalar& vs() const { return component(1); }
        scalar& vs() { return component(1); }
        
        const scalar& ws() const { return component(2); }
        scalar& ws() { return component(2); }
        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
