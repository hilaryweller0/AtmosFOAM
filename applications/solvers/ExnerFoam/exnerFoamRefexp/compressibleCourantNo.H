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

if (mesh.nInternalFaces())
{
    volScalarField CoNum = CourantNo(phi/rhof, runTime.deltaT());

    scalar maxCoNum = max(CoNum.internalField()).value();

    scalar meanCoNum = (fvc::domainIntegrate(CoNum)/gSum(mesh.V().field()))
        .value();

    scalarField sumPhi = fvc::surfaceSum(mag(phi*gUnitNormal))().internalField()
          / rho.internalField();
    scalar CoV = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    sumPhi = fvc::surfaceSum(mag(phi*(1-mag(gUnitNormal))))().internalField()
          / rho.internalField();
    scalar CoH = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    volScalarField Ndt
    (
        "Ndt",
        runTime.deltaT()*Foam::sqrt(max
        (
            -(g & fvc::grad(thetaf))/theta,
            dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
        ))
    );
    scalar maxNdt = gMax(Ndt);

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << maxCoNum << " horizontal " << CoH << " vertical " << CoV
        << " max Ndt = " << maxNdt << endl;
}

// ************************************************************************* //
