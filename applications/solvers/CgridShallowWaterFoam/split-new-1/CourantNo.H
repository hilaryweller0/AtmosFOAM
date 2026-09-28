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
    // Advective Courant number
    const scalarField CoNum
    (
        0.5*runTime.deltaTValue()/mesh.V().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()/h.primitiveField()
    );
    
    scalar CoMax = gMax(CoNum);
    scalar CoMean = gSum(CoNum);
    
    Info<< "Advective Courant number mean: " << CoMean << " max: " << CoMax << endl;

    // Gravity wave Courant number
    surfaceScalarField gwaveCoNum = dt*sqrt(linearInterpolate(magg*h))
                                    *mesh.deltaCoeffs();
    
    CoMax = gMax(gwaveCoNum.primitiveField());
    CoMean = gMax(gwaveCoNum.primitiveField());
    
    
    Info<< "Gravity wave Courant number mean: " << CoMean << " max: " << CoMax << endl;
}


// ************************************************************************* //
