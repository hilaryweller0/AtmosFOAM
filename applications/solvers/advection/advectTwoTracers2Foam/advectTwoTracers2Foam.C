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
    const bool BryanFritschSource = mesh.solutionDict().lookupOrDefault<bool>
    (
        "BryanFritschSource", true
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
            //const dimensionedScalar rhoAir(dict.lookup("rhoAir"));
        
            //rho = rho12 + rhoAir;
            fluxOld = fvc::interpolate(rho.oldTime())*phi;
            flux = fvc::interpolate(rho)*phi;
            
            if (BryanFritschSource)
            {
                S = ( 0.5*( q2-rvs + mag(q2-rvs) ) - min( ( mag(q2-rvs) - (q2-rvs) )/(mag(q2-rvs)+0.000000000001),0*rvs+1 )*min(-(q2-rvs),q1) )*transferTerm;
            }
            else
            {
                //S = (q2-JahnScale*es*T/Rv)*(q2-JahnScale*es*T/Rv) + q1*q1;
                S = (q2-JahnScale*es/(Rv*T)) - q1 + sqr((q2-JahnScale*es/(Rv*T))*(q2-JahnScale*es/(Rv*T)) + q1*q1);
                Info << "Transfer Term Min: " << min(S.internalField()).value() << " Transfer Term Max: " << max(S.internalField()).value() << endl;
            }

            // Set up the matrix without adding implicit/explicit parts
            // of advection or source terms
            fvScalarMatrix q1Eqn
            (
                //fvm::ddt(rho,q1)
                fvm::ddt(q1)
                + 0.5*fvc::div(fluxOld,q1.oldTime())
                //- 0.5*transferTerm*(q2.oldTime()-min(q1.oldTime(),rvs))
                - 0.5*timeScale*S.oldTime()
            );
            
            fvScalarMatrix q2Eqn
            (
                //fvm::ddt(rho,q2)
                fvm::ddt(q2)
                + 0.5*fvc::div(fluxOld,q2.oldTime())
                //+ 0.5*transferTerm*(q2.oldTime()-min(q1.oldTime(),rvs))
                + 0.5*timeScale*S.oldTime()
            );
            
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho)
                + 0.5*fvc::div(phi, rho.oldTime())
            );

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                q1Eqn += 0.5*fvm::div(flux,q1);
                rhoEqn += 0.5*fvm::div(phi, rho);
            }
            else
            {
                q1Eqn += 0.5*fvc::div(flux,q1);
                q2Eqn += 0.5*fvc::div(flux,q2);
                rhoEqn += 0.5*fvc::div(phi, rho);
            }
            
            // Add the damping terms, either implicit or explicit
            if (implicitSource)
            {
            //    qEqn += 0.5*fvm::Sp(rho1Damping*rho, q)
            //          + 0.5*fvm::Sp(rho2Damping*rho, q)
            //          - 0.5*rho2Damping*( rho - rhoAir );
            }
            else
            {
                //q1Eqn += -0.5*transferTerm*(q2-min(rvs,q1));
                //q2Eqn += +0.5*transferTerm*(q2-min(rvs,q1));
                q1Eqn += -0.5*timeScale*S;
                q2Eqn += +0.5*timeScale*S;
            }
            
            // Solve the matrices for the equations
            q1Eqn.solve();
            q2Eqn.solve();
            rhoEqn.solve();
            
            q2_analytic = min((rho-rhoAir)/rhoAir,q2);
            q1_analytic = (rho-rhoAir)/rhoAir - q2_analytic;
            
            Info << "Writing rho and q after rho.solve" << endl;
            
            
        }
        Info << " rho goes from " << min(rho.internalField()).value() << " to "
             << max(rho.internalField()).value() << endl;
        Info << " q1 goes from " << min(q1.internalField()).value() << " to "
             << max(q1.internalField()).value() << endl;
        Info << " Total rho in system: " << sum(rho.internalField()-1) << endl;
        Info << " rho1 fraction: " << sum(q1.internalField()*rho.internalField())/sum(rho.internalField())
         << endl;
        Info << "Transfer Term Min: " << min(S.internalField()).value() << " Transfer Term Max: " << max(S.internalField()).value() << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
