/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

{
    const scalarField sumFlux
    (
        fvc::surfaceSum(mag(flux))().primitiveField()/hf
    );
    
    const scalar CoNum =
    0.5*gMax(sumFlux/mesh.V().primitiveField())*runTime.deltaTValue();
    
    const scalar meanCoNum =
    0.5
    *(gSum(sumFlux)/gSum(mesh.V().primitiveField()))
    *runTime.deltaTValue();
    
    // Gravity wave Courant number
    surfaceScalarField gwaveCoNum = dt*sqrt(magg*hf)*mesh.deltaCoeffs();
    
    const scalar waveCoNum = gMax(gwaveCoNum.primitiveField());
    
    Info<< "Courant number mean: " << meanCoNum
    << " max: " << CoNum << endl;
    
    Info<< "Gravity wave Courant number max: " << waveCoNum
    << endl;
}


// ************************************************************************* //
