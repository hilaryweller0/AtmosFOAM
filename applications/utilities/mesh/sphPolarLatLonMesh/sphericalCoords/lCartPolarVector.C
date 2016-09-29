// The FOAM Project // File: lCartPolarVector.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   lCartPolarVector
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    stores a vector in local cartesian co-ordinates on the sphere
    and us, vs, ws (being the velocity components in the lon, lat, r directions)

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/

#include "lCartPolarVector.H"
#include "polarVector.H"
#include "polarPoint.H"
//#include "physicalConstants.H"

//using namespace Foam::physicalConstant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::lCartPolarVector::lCartPolarVector
(
    const scalar us,
    const scalar vs,
    const scalar ws 
)
{
    component(0) = us;
    component(1) = vs;
    component(2) = ws;
}        


// ************************************************************************* //
