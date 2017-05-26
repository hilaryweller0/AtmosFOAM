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
    advectTwoTracers2Foam

Description
    Solves a transport equation for two tracers rho1 and rho2 in a temperature field.
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
        
        for (int corr=0; corr < 3; corr++)
        {
            // Setup the matrix without adding implicit/explicit parts
            // of advection or source terms
            Info << "here1" << endl;
            fvScalarMatrix rho1Eqn
            (
                fvm::ddt(rho1)
              + 0.5*fvc::div(phi, rho1.oldTime())
              + 0.5*rho1damping*rho1.oldTime()
              - 0.5*rho2damping*rho2.oldTime() - 0.5*rho2damping*rho2
            );
            fvScalarMatrix rho2Eqn
            (
                fvm::ddt(rho2)
              + 0.5*fvc::div(phi, rho2.oldTime())
              + 0.5*rho2damping*rho2.oldTime()
              - 0.5*rho1damping*rho1.oldTime() - 0.5*rho1damping*rho1
            );
            fvScalarMatrix rhoEqn(fvm::ddt(rho) + 0.5*fvc::div(phi, rho.oldTime()));
            Info << "here2" << endl;

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                rho1Eqn += 0.5*fvm::div(phi, rho1);
                rho2Eqn += 0.5*fvm::div(phi, rho2);
                rhoEqn += 0.5*fvm::div(phi, rho);
            }
            else
            {
                rho1Eqn += 0.5*fvc::div(phi, rho1);
                rho2Eqn += 0.5*fvc::div(phi, rho2);
                rhoEqn += 0.5*fvc::div(phi, rho);
            }
            
            // Add the damping terms, either implicit or explicit
            if (implicitSource)
            {
                rho1Eqn += fvm::Sp(0.5*rho1damping, rho1);
                rho2Eqn += fvm::Sp(0.5*rho2damping, rho2);
            }
            else
            {
                rho1Eqn += 0.5*rho1damping*rho1;
                rho2Eqn += 0.5*rho2damping*rho2;
            }
            
            
            // Solve the matrices for the equations
            rho1Eqn.solve();
            rho2Eqn.solve();
            rhoEqn.solve();
        }
        
        Info << " rho goes from " << min(rho.internalField()).value() << " to "
             << max(rho.internalField()).value() << endl;
        Info << " rho1 goes from " << min(rho1.internalField()).value() << " to "
             << max(rho1.internalField()).value() << endl;
        Info << " rho2 goes from " << min(rho2.internalField()).value() << " to "
             << max(rho2.internalField()).value() << endl;

        Info << " Total rho in system: " << sum(rho.internalField()) << endl;
        Info << " rho1 fraction: " << sum(rho1.internalField())/sum(rho.internalField())
         << endl;
        Info << " rho2 fraction: " << sum(rho2.internalField())/sum(rho.internalField())
         << endl;
        Info << " Total fraction: " << (sum(rho1.internalField())+sum(rho2.internalField()))/sum(rho.internalField())
         << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
