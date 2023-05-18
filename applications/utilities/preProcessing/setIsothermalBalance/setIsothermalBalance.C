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
    setIsothermalBalance

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
    #include "readhRef.H"
    #include "gh.H"
    #include "createFields.H"
      
    const dictionary& itsDict = mesh.solution().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
    const label refCell = readLabel(itsDict.lookup("refCell"));
    const scalar refExner = readScalar(itsDict.lookup("refExner"));
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Iterate lnExner to convergence (should take one out iteration)
    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "lnExner iteration " << iter << endl;
        fvScalarMatrix ExnerEqn
        (
            fvm::laplacian(Cp*T, lnExner) == fvc::div(gSf)
        );
        ExnerEqn.setReference(refCell, Foam::log(refExner));
        converged = ExnerEqn.solve(Exner.name()+"Final").nIterations() == 0;
        Info << "Exner[" << refCell << "] = " << Foam::exp(lnExner[refCell]) << endl;
    }

    // Iterate Exner to convergence (slow to converge)
    Exner == exp(lnExner);
    Exner.boundaryFieldRef() = exp(lnExner.boundaryField());
    converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Exner Iteration " << iter << endl;
        
        theta == T/Exner;
        thetaf = fvc::interpolate(theta);

        fvScalarMatrix ExnerEqn
        (
            fvm::laplacian(Cp*thetaf, Exner) == fvc::div(gSf)
        );
        ExnerEqn.setReference(refCell, refExner);
        converged = ExnerEqn.solve(Exner.name()+"Final").nIterations() == 0;
        Info << "Exner[" << refCell << "] = " << Exner[refCell] << endl;
    }

    Exner.write();
    theta.write();
    volScalarField p("p", pRef*pow(Exner, 1/kappa));
    p.write();
    Exnerg = Exner - gh/(Cp*theta);
    Exnerg.write();
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

