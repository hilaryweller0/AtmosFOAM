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

\*---------------------------------------------------------------------------*/

#include "RKfield.H"
using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef RKfield<scalar, fvPatchField, volMesh> RKvolScalarField;
    typedef RKfield<scalar, fvsPatchField, surfaceMesh> RKsurfaceScalarField;

    defineTemplateTypeName(RKsurfaceScalarField);
    defineTemplateTypeName(RKvolScalarField);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::RKfield<Type, PatchField, GeoMesh>::RKfield
(
    const IOobject& io,
    const label nSteps,
    const GeometricField<Type, PatchField, GeoMesh>& f,
    const butcherTableau& BT__
)
:
    regIOobject(io),
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(nSteps),
    BT_(BT__),
    iRK_(0)
{
    for(label i = 0; i < nSteps; i++)
    {
        PtrList<GeometricField<Type, PatchField, GeoMesh> >::set(i, f);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::RKfield<Type, PatchField, GeoMesh>::RKfield
(
    const label nSteps,
    GeometricField<Type, PatchField, GeoMesh>& f,
    const butcherTableau& BT__
)
:
    regIOobject(f),
    PtrList<GeometricField<Type, PatchField, GeoMesh> >(nSteps),
    BT_(BT__),
    iRK_(0)
{
    for(label i = 0; i < nSteps; i++)
    {
        PtrList<GeometricField<Type, PatchField, GeoMesh> >::set(i, f);
    }
}


// ************************************************************************* //
