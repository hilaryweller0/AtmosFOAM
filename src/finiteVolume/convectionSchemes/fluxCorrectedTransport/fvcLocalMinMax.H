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
    Functions to set the volume field based on the min and max of all neighbours

SourceFiles
    fvcLocalMinMax.C

\*---------------------------------------------------------------------------*/


#ifndef fvcLocalMinMax_H
#define fvcLocalMinMax_H

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
    void localMin
    (
        Field<scalar>&,
        const surfaceScalarField&
    );

    tmp<volScalarField>
    localMin
    (
        const surfaceScalarField&
    );

    tmp<volScalarField>
    localMin
    (
        const tmp<surfaceScalarField>&
    );

    void localMax
    (
        Field<scalar>&,
        const surfaceScalarField&
    );

    tmp<volScalarField>
    localMax
    (
        const surfaceScalarField&
    );

    tmp<volScalarField>
    localMax
    (
        const tmp<surfaceScalarField>&
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
