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
#include "localMinMin.H"
#include "fvcLocalMinMax.H"
#include "upwind.H"
#include "downwind.H"
//#include "uncorrectedSnGrad.H"
//#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fluxLimit
(
    volScalarField& Td,
    const surfaceScalarField& fluxCorr,
    const dimensionedScalar& dt,
    const int nIter
)
{
    const fvMesh& mesh = fluxCorr.mesh();
    const volScalarField& T0 = Td.oldTime();

    // Schemes needed
    Foam::localMax<scalar> maxInterp(mesh);
    Foam::localMinMin<scalar> minInterp(mesh);

    // Where can we use T from the old time?

    // Local extrema
    volScalarField Tmin = min(Td, T0);
    volScalarField Tmax = max(Td, T0);
    surfaceScalarField Tfmax = maxInterp.interpolate(Tmax);
    surfaceScalarField Tfmin = minInterp.interpolate(Tmin);
    Tmax = fvc::localMax(Tfmax);
    Tmin = fvc::localMin(Tfmin);

    fluxLimit(Td, fluxCorr, Tmin, Tmax, dt, nIter);
}


template<class Type>
void fluxLimit
(
    volScalarField& Td,
    const surfaceScalarField& fluxCorr,
    const Type& Tmin,
    const Type& Tmax,
    const dimensionedScalar& dt,
    const int nIter
)
{
    // Gradient schemes to check for diffusive correction special cases
    //const fvMesh& mesh = Td.mesh();
    //fv::uncorrectedSnGrad<scalar> snGrad(mesh);
    //fv::gaussGrad<scalar> gGrad(mesh);

    surfaceScalarField limitedCorr = fluxCorr;
    for(int iter = 0; iter < nIter; iter++)
    {
        if (iter > 0)
        {
            limitedCorr = fluxCorr - limitedCorr;
        }
        // Amount each cell can rise or fall by
        volScalarField Qp = Tmax - Td;
        volScalarField Qm = Td - Tmin;

        // Check for diffusive (rather than anti-diffusive) corrections
        // (not used as it makes the results noisy)
        /*surfaceScalarField gradnT = snGrad.snGrad(Td);
        volVectorField gradT = gGrad.grad(Td);
        forAll(limitedCorr, faceI)
        {
            scalar& A = limitedCorr[faceI];
            label own = mesh.owner()[faceI];
            label nei = mesh.neighbour()[faceI];
            scalar gradTo = gradT[own] & mesh.Sf()[faceI];
            scalar gradTn = gradT[nei] & mesh.Sf()[faceI];

            if (A*gradnT[faceI] < 0 && (A*gradTo < 0 || A*gradTn < 0))
            {
                //A = 0;
            }
        }*/

        fluxLimitFromQ(limitedCorr, Qp, Qm, dt);

        Td -= dt*fvc::div(limitedCorr);
    }
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
        if (mag(Pp[cellI]) < SMALL) Rp[cellI] = 0;
        if (mag(Pm[cellI]) < SMALL) Rm[cellI] = 0;
    }

    // Up and downwind interpolation schemes based on mesh orientation
    upwind<scalar> up(mesh, mesh.magSf());
    downwind<scalar> down(mesh, mesh.magSf());

    // Limit the fluxes
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
