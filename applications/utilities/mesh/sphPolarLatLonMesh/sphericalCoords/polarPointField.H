// The FOAM Project // File: polarPointField.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   polarPointField
   \\  /           | Family: polarMesh
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
TYPEDEF
    polarPointField

DESCRIPTION
    polarPointField is a Field of polarPoint s.

AUTHOR
    Hilary Spencer.

*/
// ------------------------------------------------------------------------- //

#ifndef polarPointField_H
#define polarPointField_H

#include "polarPoint.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef Field<polarPoint> polarPointField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
