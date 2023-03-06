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

Application
    setThetaExnerFromT

Description
    Find discretely balanced theta and Exner profiles given temperature

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "ExnerTheta.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermo.H"
    #include "createFields.H"
      
    const dictionary& itsDict = mesh.solutionDict().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
    const label refCell = readLabel(itsDict.lookup("refCell"));
    const scalar refExner = readScalar(itsDict.lookup("refExner"));
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Iteration " << iter << endl;
        
        theta == T/ExnerInit;
        thetaf = fvc::interpolate(theta);

        fvScalarMatrix ExnerEqn
        (
            fvm::laplacian(Cp*thetaf, ExnerInit) == fvc::div(gSf)
        );
        ExnerEqn.setReference(refCell, refExner);
        converged = ExnerEqn.solve(Exner.name()).nIterations() == 0;
        
        Info << "Exner[" << refCell << "] = " << ExnerInit[refCell] << endl;
        Info << "T[" << refCell << "] = " << T[refCell] << endl;
        Info << "theta[" << refCell << "] = " << theta[refCell] << endl;
    }
    theta.write();
    
    Exner == ExnerInit;
    Exner.boundaryFieldRef() = ExnerInit.boundaryField();
    Exner.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

