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
    Applies Zalesak (1979) flux limiter to a flux

SourceFiles
    fvcFluxLimit.C

\*---------------------------------------------------------------------------*/


#ifndef fvcFluxLimit_H
#define fvcFluxLimit_H

#include "fvcSurfaceIntegrateInOut.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    void fluxLimit
    (
        surfaceScalarField& fluxCorr,
        const volScalarField& Td,
        const volScalarField& Told,
        const dimensionedScalar& dt
    );

    void fluxLimit
    (
        surfaceScalarField& fluxCorr,
        const volScalarField& Td,
        const dimensionedScalar& dt
    );

    void fluxLimit
    (
        surfaceScalarField& fluxCorr,
        const volScalarField& Td,
        const dimensionedScalar& Tmin,
        const dimensionedScalar& Tmax,
        const dimensionedScalar& dt
    );

    void fluxLimit
    (
        surfaceScalarField& fluxCorr,
        const volScalarField& Td,
        const volScalarField& Tmin,
        const volScalarField& Tmax,
        const dimensionedScalar& dt
    );

    void fluxLimitFromQ
    (
        surfaceScalarField& fluxCorr,
        const volScalarField& Qp,
        const volScalarField& Qm,
        const dimensionedScalar& dt
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
