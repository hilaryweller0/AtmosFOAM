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
    advectTwoTracersFoam

Description
    Solves a transport equation for two tracers T1 and T2.
    CODE IS WORK IN PROGRESS

\*---------------------------------------------------------------------------*/

#include "gdal_priv.h"
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

    const bool implicitAdvection = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitAdvection", false
    );
    const bool implicitSource = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitSource", false
    );

    #include "calculateSource.H"

    Info<< "\nCalculating advection\n" << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
        if (!implicitAdvection && CoNum > 1.0)
        {
            FatalErrorInFunction << "Max Courant number > 1"
                << exit(FatalError);
        }

        for (int corr=0; corr < 3; corr++)
        {
            // Setup the matrix without adding implicit/explicit parts
            // of advection or source terms
            fvScalarMatrix T1Eqn
            (
                fvm::ddt(T1)
              + 0.5*fvc::div(phi, T1.oldTime())
              + 0.5*T1damping*T1.oldTime()
              - 0.5*T2damping*T2.oldTime() - 0.5*T2damping*T2
            );
            fvScalarMatrix T2Eqn
            (
                fvm::ddt(T2)
              + 0.5*fvc::div(phi, T2.oldTime())
              + 0.5*T2damping*T2.oldTime()
              - 0.5*T1damping*T1.oldTime() - 0.5*T1damping*T1
            );
            fvScalarMatrix TEqn(fvm::ddt(T) + 0.5*fvc::div(phi, T.oldTime()));

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                T1Eqn += 0.5*fvm::div(phi, T1);
                T2Eqn += 0.5*fvm::div(phi, T2);
                TEqn += 0.5*fvm::div(phi, T);
            }
            else
            {
                T1Eqn += 0.5*fvc::div(phi, T1);
                T2Eqn += 0.5*fvc::div(phi, T2);
                TEqn += 0.5*fvc::div(phi, T);
            }
            
            // Add the damping terms, either implicit or explicit
            if (implicitSource)
            {
                T1Eqn += fvm::Sp(0.5*T1damping, T1);
                T2Eqn += fvm::Sp(0.5*T2damping, T2);
            }
            else
            {
                T1Eqn += 0.5*T1damping*T1;
                T2Eqn += 0.5*T2damping*T2;
            }
            
            // Solve the matrices for the equations
            T1Eqn.solve();
            T2Eqn.solve();
            TEqn.solve();
        }
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        Info << " T1 goes from " << min(T1.internalField()).value() << " to "
             << max(T1.internalField()).value() << endl;
        Info << " T2 goes from " << min(T2.internalField()).value() << " to "
             << max(T2.internalField()).value() << endl;

        runTime.write();
    }
    
    Info << " Total T in system: " << sum(T.internalField()) << endl;
    Info << " T1 fraction: " << sum(T1.internalField())/sum(T.internalField())
         << endl;
    Info << " T2 fraction: " << sum(T2.internalField())/sum(T.internalField())
         << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
