/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

Class
    Foam::RKfield

Description
    Holds a GeometricField at a number of intermediate steps for use by
    Runge-Kutta time integration with a doubleButcherTableau

SourceFiles
    RKfieldI.H
    RKfield.C
    RKfieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef RKfield_H
#define RKfield_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class RKfield Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class RKfield
:
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        RKfield(const RKfield<Type, PatchField, GeoMesh>&);

        //- Disallow default bitwise assignment
        void operator=(const RKfield<Type, PatchField, GeoMesh>&);


public:

    // Constructors

        //- Construct from size and one field
        RKfield
        (
            const label nSteps,
            const GeometricField<Type, PatchField, GeoMesh>& f
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "RKfield.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
