/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    continuityErrs

Description
    Calculates and prints the continuity errors.

\*---------------------------------------------------------------------------*/

{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    volScalarField rhoDiff("rhoDiff", rho - Psi*Exner);

    scalar sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rhoDiff))/totalMass
    ).value();

    scalar globalContErr =
    (
        fvc::domainIntegrate(rhoDiff)/totalMass
    ).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
    
    Info << "Maximum rho difference fraction = "
         //<< max(mag(rhoDiff.internalField()/rho.internalField())) << endl;
         << max(mag(rhoDiff/rho)).value() << endl;

//    // Finding location of maximum U and Co there
//    scalar Umax = 0;
//    label cellMax = -1;
//    forAll(U, cellI)
//    {
//        if (mag(U[cellI]) > Umax)
//        {
//            Umax = mag(U[cellI]);
//            cellMax = cellI;
//        }
//    }
//    Info << "Max U " << Umax << " at " << cellMax << " rho = " << rho[cellMax]
//         << " V = " << mesh.V()[cellMax] << "\nFaces, areas and fluxes:\n";
//    for(label i = 0; i < mesh.cells()[cellMax].size(); i++)
//    {
//        label faceI = mesh.cells()[cellMax][i];
//        Info << faceI << "  " << mag(mesh.Sf()[faceI]) << "  " << phi[faceI] << nl;
//    }
}

// ************************************************************************* //
