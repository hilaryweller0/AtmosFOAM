/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifndef setInternalValues_H
#define setInternalValues_H

#include "GeometricFields.H"

#define TEMPLATE \
    template<class vector, template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

// Set values of f1 to those of f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const GeometricField<Type, PatchField, GeoMesh>& f2
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorIn("setInternalValues") << "Cannot set values of field "
            << f1.name() << " with size " << f1.size() << " to values of field"
            << f2.name() << " with size " << f2.size() << abort(FatalError);
    }
    
    f1.ref().Field<Type>::operator=(f2.internalField());
}

// Set values of f1 to those of f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalAndBoundaryValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const GeometricField<Type, PatchField, GeoMesh>& f2
)
{
    setInternalValues(f1, f2);

    // Also set all the boundary values
    for(label patchI = 0; patchI < f1.boundaryField().size(); patchI++)
    {
        f1.boundaryFieldRef()[patchI] == f2.boundaryField()[patchI];
    }
}


// Set values of f1 to those of f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalAndBoundaryValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tf2
)
{
    const GeometricField<Type, PatchField, GeoMesh>& f2 = tf2();
    setInternalValues(f1, f2);

    // Also set all the boundary values
    for(label patchI = 0; patchI < f1.boundaryField().size(); patchI++)
    {
        f1.boundaryFieldRef()[patchI] == f2.boundaryField()[patchI];
    }
    tf2.clear();
}


// Set values of f1 to  f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const Field<Type>& f2
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorIn("setInternalValues") << "Cannot set values of field "
            << f1.name() << " with size " << f1.size() << " to values of field"
            << " of size " << f2.size() << abort(FatalError);
    }
    
    f1.ref().Field<Type>::operator=(f2);
}


// Set values of f1 to  f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const DimensionedField<Type, GeoMesh>& f2
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorIn("setInternalValues") << "Cannot set values of field "
            << f1.name() << " with size " << f1.size() << " to values of field"
            << " of size " << f2.size() << abort(FatalError);
    }
    
    f1.ref().Field<Type>::operator=(f2);
}

// Set values of f1 to  f2
template<class Type, template<class> class PatchField, class GeoMesh>
inline void setInternalValues
(
    GeometricField<Type, PatchField, GeoMesh>& f1,
    const tmp<DimensionedField<Type, GeoMesh> >& tf2
)
{
    const DimensionedField<Type, GeoMesh>& f2 = tf2();
    setInternalValues(f1, f2);
    tf2.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#endif
// ************************************************************************* //
