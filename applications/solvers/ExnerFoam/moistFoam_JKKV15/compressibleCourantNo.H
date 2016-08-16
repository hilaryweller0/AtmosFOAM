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
    Calculates and outputs the mean and maximum Courant Numbers and Ndt

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar maxNdt = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(U))().internalField()
      / rho.internalField()
    );

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    volScalarField Ndt
    (
        "Ndt",
        dt*Foam::sqrt(mag((g & fvc::grad(thetaRho))/thetaRho))
    );
    maxNdt = gMax(Ndt);
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << " max Ndt = " << maxNdt << endl;

// ************************************************************************* //
