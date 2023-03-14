/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "HodgeOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HodgeOps::HodgeOps(const fvMesh& mesh__)
:
    mesh_(orthogonalBoundaries(mesh__)),
    delta_(mesh_.delta()),
    magd_(mag(delta_)),
    Hdiag_((mesh_.Sf() & delta_)/sqr(magd_))
{
    surfaceScalarField dc = mesh_.deltaCoeffs();
    const_cast<surfaceScalarField&>(mesh_.nonOrthDeltaCoeffs()) = dc;
}

// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::HodgeOps::orthogonalBoundaries(const fvMesh& mesh__)
{
    Foam::Info << "Making the boundary faces orthogonal" << Foam::endl;

    forAll(mesh__.boundary(), patchi)
    {
        const fvPatch& pat = mesh__.boundary()[patchi];
        // only move centres of un-coupled boundary faces
        if (!pat.coupled())
        {
            // loop through and move all of the faces
            for(label facei = pat.start(); facei <pat.start()+pat.size();facei++)
            {
                label own = mesh__.faceOwner()[facei];
                const vector& C = mesh__.C()[own];
                const vector& Cf = mesh__.faceCentres()[facei];
                const vector Sfhat = mesh__.faceAreas()[facei]
                                      /mag(mesh__.faceAreas()[facei]);
                const vector fCnew = C + ((Cf - C) & Sfhat)*Sfhat;

                const_cast<vectorField&>(mesh__.faceCentres())[facei] = fCnew;
            }
        }
    }

    return mesh__;
}

// ************************************************************************* //

