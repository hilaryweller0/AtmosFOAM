/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the maximum Courant Number
    and the maximum gravity wave Courant number

\*---------------------------------------------------------------------------*/

scalar maxCoG = 0.0;
scalar maxCo = 0.0;

if (mesh.nInternalFaces())
{
    maxCoG = max
    (
        runTime.deltaTValue()*mesh.deltaCoeffs()*Foam::sqrt(g*hc)
    ).value();

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(volFlux))().internalField()
    );

    maxCo = 0.5*gMax(sumPhi/(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Max gravirty wave Courant number: " << maxCoG
    << " advective Cournat number maximum: " << maxCo << endl;

// ************************************************************************* //
