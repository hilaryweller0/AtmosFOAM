// The FOAM Project // File: polarPoint.C
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

#include "polarPoint.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::polarPoint::degToRad = Foam::constant::mathematical::pi/180.;

Foam::scalar Foam::polarPoint::radToDeg = 180./Foam::constant::mathematical::pi;

// ************************************************************************* //
