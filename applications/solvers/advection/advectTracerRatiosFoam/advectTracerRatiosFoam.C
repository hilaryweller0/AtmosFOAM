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
            const dimensionedScalar unitDamping("unitDamping", dimensionSet(0,0,0,0,0), 1);
            q2 = -q+1;
            fluxOld = fvc::interpolate(T.oldTime())*phi;
            flux = fvc::interpolate(T)*phi;
            // Setup the matrix without adding implicit/explicit parts
            // of advection or source terms
            fvScalarMatrix qEqn
            (
                //fvm::ddt(q)
              fvm::ddt(T,q)
              + 0.5*fvc::div(fluxOld,q.oldTime())
              //+ 0.5*( U & fvc::grad(q.oldTime()) )
              + 0.5*T1damping*q.oldTime()
              - 0.5*T2damping*(1-q.oldTime())
            );
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + 0.5*fvc::div(phi, T.oldTime())
            );

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                qEqn += 0.5*fvm::div(flux,q);
                TEqn += 0.5*fvm::div(phi, T);
            }
            else
            {
                qEqn += 0.5*fvc::div(flux,q);
                //qEqn += 0.5*( U & fvc::grad(q) );
                TEqn += 0.5*fvc::div(phi, T);
            }
            
            // Add the damping terms, either implicit or explicit
            if (implicitSource)
            {
                //qEqn += 0.5*fvm::Sp(T1damping, q);
                //qEqn += -0.5*fvm::Sp(T2damping, q2);
                //qEqn += 0.5*fvm::Sp(unitDamping, q*T1damping);
                //qEqn += -0.5*fvm::Sp(unitDamping, q2*T2damping);
            }
            else
            {
                qEqn += 0.5*T1damping*q;
                qEqn += -0.5*T2damping*q2;
            }
            
            // Solve the matrices for the equations
            TEqn.solve();
            Info << "Writing T and q after T.solve" << endl;
            T.write();
            q.write();
            qEqn.solve();
            
        }
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        Info << " q goes from " << min(q.internalField()).value() << " to "
             << max(q.internalField()).value() << endl;
        Info << " Total T in system: " << sum(T.internalField()) << endl;
        Info << " T1 fraction: " << sum(q.internalField()*T.internalField())/sum(T.internalField())
         << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
