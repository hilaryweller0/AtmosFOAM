/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    advectTracerRatiosFoam

Description
    Solves a transport equation for two tracers T1 and T2.
    CODE IS WORK IN PROGRESS

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "velocityField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        
        for (int corr=0; corr < 3; corr++)
        {
            fvScalarMatrix rhoaEqn
            (
                fvm::ddt(rho_a)
                + 0.5*fvc::div(phi, rho_a.oldTime())
                + 0.5*fvc::div(phi, rho_a)
            );
            
            fvScalarMatrix rhobEqn
            (
                fvm::ddt(rho_b)
                + 0.5*fvc::div(phi, rho_b.oldTime())
                + 0.5*fvc::div(phi, rho_b)
            );
            
            // Solve the matrices for the equations
            rhoaEqn.solve();
            rhobEqn.solve();
        }
        Info << " rho_a goes from " << min(rho_a.internalField()).value() << " to "
             << max(rho_a.internalField()).value() << endl;
        Info << " rho_b goes from " << min(rho_b.internalField()).value() << " to "
             << max(rho_b.internalField()).value() << endl;
        runTime.write();
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
