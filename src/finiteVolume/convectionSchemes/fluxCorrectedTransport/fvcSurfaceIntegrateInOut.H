/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InNamespace
    Foam::fvc

Description
    Surface integrate surfaceField creating a volField.
    Separate sums for fluxes into an out of a cell

SourceFiles
    fvcSurfaceIntegrateInOut.C

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceIntegrateInOut_H
#define fvcSurfaceIntegrateInOut_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    void surfaceIntegrateIn
    (
        Field<Type>&,
        const Foam::SurfaceField<Type>&
    );

    template<class Type>
    tmp<Foam::VolField<Type>>
    surfaceIntegrateIn
    (
        const Foam::SurfaceField<Type>&
    );

    template<class Type>
    tmp<Foam::VolField<Type>>
    surfaceIntegrateIn
    (
        const tmp<Foam::SurfaceField<Type>>&
    );

    template<class Type>
    void surfaceIntegrateOut
    (
        Field<Type>&,
        const Foam::SurfaceField<Type>&
    );

    template<class Type>
    tmp<Foam::VolField<Type>>
    surfaceIntegrateOut
    (
        const Foam::SurfaceField<Type>&
    );

    template<class Type>
    tmp<Foam::VolField<Type>>
    surfaceIntegrateOut
    (
        const tmp<Foam::SurfaceField<Type>>&
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "fvcSurfaceIntegrateInOut.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
