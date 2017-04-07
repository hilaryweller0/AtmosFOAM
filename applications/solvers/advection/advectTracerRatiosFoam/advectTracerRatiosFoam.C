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
        //if (!implicitAdvection && CoNum > 1.0)
        //{
        //    FatalErrorInFunction << "Max Courant number > 1"
        //        << exit(FatalError);
        //}

        for (int corr=0; corr < 3; corr++)
        {
            const dimensionedScalar rhoAir(dict.lookup("rhoAir"));
        
            //rho = rho12 + rhoAir;
            fluxOld = fvc::interpolate(rho.oldTime())*phi;
            flux = fvc::interpolate(rho)*phi;
            
            // Set up the matrix without adding implicit/explicit parts
            // of advection or source terms
            fvScalarMatrix qEqn
            (
                fvm::ddt(rho,q)
                + 0.5*fvc::div(fluxOld,q.oldTime())
                + 0.5*rho1Damping*q.oldTime()*rho.oldTime()
                - 0.5*rho2Damping*( (1-q.oldTime())*rho.oldTime() - rhoAir )
                //+ 0.5*rho1Damping*0.5*( q.oldTime()*rho.oldTime() + mag(q.oldTime()*rho.oldTime()) )
                //- 0.5*rho2Damping*0.5*( (1 - q.oldTime())*rho.oldTime() - rhoAir + mag( (1 - q.oldTime())*rho.oldTime() - rhoAir ) )
            );
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho)
                + 0.5*fvc::div(phi, rho.oldTime())
            );

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                qEqn += 0.5*fvm::div(flux,q);
                rhoEqn += 0.5*fvm::div(phi, rho);
            }
            else
            {
                qEqn += 0.5*fvc::div(flux,q);
                rhoEqn += 0.5*fvc::div(phi, rho);
            }
            
            // Add the damping terms, either implicit or explicit
            if (implicitSource)
            {
                qEqn += 0.5*fvm::Sp(rho1Damping*rho, q)
                      + 0.5*fvm::Sp(rho2Damping*rho, q)
                      - 0.5*rho2Damping*( rho - rhoAir );
            }
            else
            {
                //qEqn += 0.5*rho1Damping*0.5*( q*rho + mag(q*rho) )
                //      - 0.5*rho2Damping*0.5*( (1 - q)*rho - rhoAir + mag( (1 - q)*rho - rhoAir ) );
                
                
                //qEqn += 0.5*rho1Damping*q*rho;
                qEqn += - 0.5*rho2Damping*((1 - q)*rho - rhoAir);
                
                qEqn += 0.5*fvm::Sp(rho1Damping*rho, q);
                //qEqn += 0.5*fvm::Sp(rho2Damping*rho, q)
                //      - 0.5*rho2Damping*( rho - rhoAir );
            }
            
            // Solve the matrices for the equations
            
            qEqn.solve();
            rhoEqn.solve();
            
            Info << "Writing rho and q after rho.solve" << endl;
            
            
        }
        Info << " rho goes from " << min(rho.internalField()).value() << " to "
             << max(rho.internalField()).value() << endl;
        Info << " q goes from " << min(q.internalField()).value() << " to "
             << max(q.internalField()).value() << endl;
        Info << " Total rho in system: " << sum(rho.internalField()-1) << endl;
        Info << " rho1 fraction: " << sum(q.internalField()*rho.internalField())/sum(rho.internalField())
         << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
