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
    Calculates and outputs maximum Courant Numbers for each partition

\*---------------------------------------------------------------------------*/

scalarList maxCoNums(nParts, scalar(0));

if (mesh.nInternalFaces())
{
    for(label ip = 0; ip < nParts; ip++)
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(flux[ip]))().internalField()
        );

        maxCoNums[ip] = 0.5*gMax(sumPhi/(mesh.V().field()))
                    *runTime.deltaTValue();
    }
}

Info<< "Courant partition maximums: " << maxCoNums << endl;

// ************************************************************************* //
