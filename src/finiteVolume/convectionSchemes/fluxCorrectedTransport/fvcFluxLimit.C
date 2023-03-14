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

\*---------------------------------------------------------------------------*/

#include "fvcFluxLimit.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "localMax.H"
#include "localMin.H"
#include "fvcLocalMinMax.H"
#include "upwind.H"
#include "downwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fluxLimit
(
    surfaceScalarField& fluxCorr,
    const volScalarField& Td,
    const volScalarField& Told,
    const dimensionedScalar& dt
)
{
    const fvMesh& mesh = fluxCorr.mesh();
    
    // Local extrema
    volScalarField Tmin = min(Td, Told);
    volScalarField Tmax = max(Td, Told);
    Foam::localMax<scalar> maxInterp(mesh);
    Foam::localMin<scalar> minInterp(mesh);
    surfaceScalarField Tfmax = maxInterp.interpolate(Tmax);
    surfaceScalarField Tfmin = minInterp.interpolate(Tmin);
    Tmax = fvc::localMax(Tfmax);
    Tmin = fvc::localMin(Tfmin);

    // Amount each cell can rise or fall by
    volScalarField Qp = Tmax - Td;
    volScalarField Qm = Td - Tmin;

    fluxLimitFromQ(fluxCorr, Qp, Qm, dt);
}


void fluxLimit
(
    surfaceScalarField& fluxCorr,
    const volScalarField& Td,
    const dimensionedScalar& dt
)
{
    const fvMesh& mesh = fluxCorr.mesh();
    
    // Local extrema
    Foam::localMax<scalar> maxInterp(mesh);
    Foam::localMin<scalar> minInterp(mesh);
    surfaceScalarField Tfmax = maxInterp.interpolate(Td);
    surfaceScalarField Tfmin = minInterp.interpolate(Td);
    volScalarField Tmax = fvc::localMax(Tfmax);
    volScalarField Tmin = fvc::localMin(Tfmin);

    // Amount each cell can rise or fall by
    volScalarField Qp = Tmax - Td;
    volScalarField Qm = Td - Tmin;

    fluxLimitFromQ(fluxCorr, Qp, Qm, dt);
}


void fluxLimit
(
    surfaceScalarField& fluxCorr,
    const volScalarField& Td,
    const dimensionedScalar& Tmin,
    const dimensionedScalar& Tmax,
    const dimensionedScalar& dt
)
{
    // Amount each cell can rise or fall by
    volScalarField Qp = Tmax - Td;
    volScalarField Qm = Td - Tmin;

    fluxLimitFromQ(fluxCorr, Qp, Qm, dt);
}


void fluxLimitFromQ
(
    surfaceScalarField& fluxCorr,
    const volScalarField& Qp,
    const volScalarField& Qm,
    const dimensionedScalar& dt
)
{
    const fvMesh& mesh = fluxCorr.mesh();
    
    // Fluxes into and out of each cell
    volScalarField Pp = dt*fvc::surfaceIntegrateIn(fluxCorr);
    volScalarField Pm = dt*fvc::surfaceIntegrateOut(fluxCorr);
    
    // Ratios
    dimensionedScalar Psmall("Psmall", Pp.dimensions(), SMALL);
    volScalarField Rp = min(1., Qp/max(Pp, Psmall));
    volScalarField Rm = min(1., Qm/max(Pm, Psmall));
    forAll(Rp, cellI)
    {
        if (mag(Pp[cellI]) < VSMALL) Rp[cellI] = 0;
        if (mag(Pm[cellI]) < VSMALL) Rm[cellI] = 0;
    }
    
    // Limit the fluxes
    upwind<scalar> up(mesh, mesh.magSf());
    downwind<scalar> down(mesh, mesh.magSf());
    surfaceScalarField RpOwn = up.interpolate(Rp);
    surfaceScalarField RpNei = down.interpolate(Rp);
    surfaceScalarField RmOwn = up.interpolate(Rm);
    surfaceScalarField RmNei = down.interpolate(Rm);
    
    fluxCorr *= 0.5*(1+sign(fluxCorr))*min(RpNei, RmOwn)
              + 0.5*(1-sign(fluxCorr))*min(RpOwn, RmNei);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
