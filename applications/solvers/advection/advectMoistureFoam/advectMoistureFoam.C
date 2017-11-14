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
    Advection of moisture in 2 ways.
    Diagnostic: Advect total moisture and diagnose liquid and vapour based on local temperature.
    Prognostic: Advect liquid and vapour with phase changes between the two.

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
    #include "readDicts.H"
    #include "createFields.H"

    bool implicitSource = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitSource", false
    );
    const bool implicitConservation = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitConservation", false
    );
    double offCentre = mesh.solutionDict().lookupOrDefault<double>
    (
        "CrankNicolsonAlpha", 0.5
    );
    const int numberOfCorrectors = mesh.solutionDict().lookupOrDefault<int>
    (
        "numberOfCorrectors", 3
    );
    const int opSplit = mesh.solutionDict().lookupOrDefault<int>
    (
        "operatorSplitting", 0
    );
    const double opSplitOffCentre = mesh.solutionDict().lookupOrDefault<double>
    (
        "opSplitOffCenterCoefficient", 0
    );
    const double timeScale = mesh.solutionDict().lookupOrDefault<double>
    (
        "timeScale", 1.
    );
    if (opSplit)
    {
        implicitSource = false;
        offCentre = 0.;
    }

    #include "writeDiagnosticsInit.H"

    Info<< "\nCalculating advection\n" << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
        
        rl_temp = rl;
        rv_temp = rv;
        
        for (int corr=0; corr < numberOfCorrectors; corr++)
        {
            rl = rl.oldTime();
            rv = rv.oldTime();
        
            #include "calculateSource.H"
            //Operator Splitting, pre-advection step.
            if (opSplit)
            {
                rl -= dt*(1-opSplitOffCentre)*(Sl*rl_temp - Sv*(rv_temp-rvs));
                rv += dt*(1-opSplitOffCentre)*(Sl*rl_temp - Sv*(rv_temp-rvs));
                rl_temp = rl;
                rv_temp = rv;
            }
            
            
            //Prognostic Method.
            rl -= dt*
            (
                (1-offCentre)*fvc::div(phi,rl.oldTime())
              +     offCentre*fvc::div(phi,rl_temp)
              + (1-opSplit)*
                (
                    (1-offCentre)*Sl*rl.oldTime()
                  - (1-offCentre)*Sv*(rv.oldTime()-rvs)
                  -     offCentre*Sv*(rv_temp-rvs)
                )
            );
            if (implicitSource) rl =  rl/(dt*offCentre*Sl+1);
            else                rl -= (1-opSplit)*dt*offCentre*Sl*rl_temp;
            
            rv -= dt* 
            (
                (1-offCentre)*fvc::div(phi,rv.oldTime())
              +     offCentre*fvc::div(phi,rv_temp)
              + (1-opSplit)*
                (
                  - (1-offCentre)*Sl*rl.oldTime()
                  + (1-offCentre)*Sv*(rv.oldTime()-rvs)
                  -     offCentre*Sv*rvs
                )
            );
            if (implicitSource)
            {
                rv += dt*offCentre*Sl*rl;
                rv =  rv/(dt*offCentre*Sv+1);
                if (implicitConservation)
                    rl -= dt*offCentre*Sv*(rv_temp-rv);
            }
            else                
            {
                rv += (1-opSplit)*dt*offCentre*(Sl*rl_temp - Sv*rv_temp);
            }
            
            
            //Operator Splitting, post-advection step.
            if (opSplit)
            {
                #include "calculateSource.H"
                rl_temp = rl;
                rv_temp = rv;
                rl -= dt*opSplitOffCentre*(Sl*rl_temp - Sv*(rv_temp-rvs));
                rv += dt*opSplitOffCentre*(Sl*rl_temp - Sv*(rv_temp-rvs));
            }
                        
            
            rl_temp = rl;
            rv_temp = rv;
            
            
            //Diagnostic Method
            fvScalarMatrix rtEqn
            (
                fvm::ddt(rt)
                + (1-offCentre)*fvc::div(phi, rt.oldTime())
                +     offCentre*fvc::div(phi, rt)
            );
            rtEqn.solve();
            
            rv_diag = min(rt,rvs);
            rl_diag = rt - rv_diag;
        
            #include "analyticSolution.H" 
        }
        
        #include "writeDiagnostics.H"
        
        Info << " rt goes from " 
             << min(rt.internalField()).value() << " to "
             << max(rt.internalField()).value() << endl;
        Info << " rl goes from " 
             << min(rl.internalField()).value() << " to "
             << max(rl.internalField()).value() << endl;
        Info << " Total rt in system: " 
             << sum(rt.internalField()-1) << endl;
        Info << " Total Moisture: " 
             << sum(rl.internalField())+sum(rv.internalField()) << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
