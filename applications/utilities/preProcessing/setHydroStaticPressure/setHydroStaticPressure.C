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
    setHydroStaticPressure

Description
    Find discretely hydrostacally balanced p_rgh in balance with given theta

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "ExnerTheta.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "readThermoProperties.H"
    #include "createFields.H"
      
    const dictionary& itsDict = mesh.solution().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
    const int BCiters  = itsDict.lookupOrDefault<int>("BCiters", 10);
    const scalar BCtol = itsDict.lookupOrDefault<scalar>("BCtol", SMALL);
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Find top and bottom boundaries
    label groundBC = -1;
    label topBC = -1;
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (mesh.boundaryMesh()[patchi].name() == "top")
        {
            topBC = patchi;
        }
        else if (mesh.boundaryMesh()[patchi].name() == "ground")
        {
            groundBC = patchi;
        }
    }
    if (groundBC == -1 || topBC == -1)
    {
        FatalErrorIn("setHydroStaticPressure")
            << " no boundary called top or ground. The boundaries are called "
            << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    // Iterate until boundaries are as required and in hydrostatic balance
    bool innerConverged = false;
    bool outerConverged = false;
    scalar topBCval = p_rgh.boundaryField()[topBC][0];
    for
    (
        label iter = 0;
        iter < maxIters && !(innerConverged && outerConverged);
        iter++
    )
    {
        Info << "Outer iteration " << iter << endl;
        
        // Inner iterations with fixed top boundary
        innerConverged = false;
        for(label BCiter = 0; BCiter < BCiters && !innerConverged; BCiter++)
        {
            p == p_rgh + rho*gh;
            thermo.T() == pow(p/pRef, kappa)*theta;
            thermo.he() == thermo.he(p,thermo.T());
            thermo.correct();
            rho == thermo.rho();
            phi == -gzSf*fvc::snGrad(rho);

            fvScalarMatrix p_rghEqn
            (
                fvc::div(phi)
              - fvm::laplacian(p_rgh)
            );
            innerConverged = p_rghEqn.solve(p_rgh.name())
                    .nIterations() == 0;

        }
        scalar maxGroundP_rgh = max(p_rgh.boundaryField()[groundBC]);
        outerConverged = (mag(pRef.value()-maxGroundP_rgh)< BCtol);
        
        // modify the top boundary condition
        Info << "p_rgh ground value = " << maxGroundP_rgh
             << "  ground value minus pRef = "
             << pRef.value()-maxGroundP_rgh
             << " p_rgh top BC going from  = " << topBCval << " to ";

        topBCval *= pRef.value()/maxGroundP_rgh;
        topBCval = min(max(topBCval, scalar(0)), pRef.value());
        Info << topBCval << endl;
        p_rgh.boundaryFieldRef()[topBC] == topBCval;
    }
    
    p_rgh.write();
    p.write();
    thermo.T().write();
    rho.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

