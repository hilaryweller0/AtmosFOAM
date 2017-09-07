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
    #include "createFields.H"

    const bool implicitAdvection = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitAdvection", false
    );
    const bool implicitSource = mesh.solutionDict().lookupOrDefault<bool>
    (
        "implicitSource", false
    );
    const double a = mesh.solutionDict().lookupOrDefault<double>
    (
        "CrankNicolsonAlpha", 0.5
    );

    #include "calculateSource.H"
    #include "writeDiagnosticsInit.H"

    Info<< "\nCalculating advection\n" << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
        
        rl_temp = rl;
        rv_temp = rv;
        
        for (int corr=0; corr < 3; corr++)
        {
            //const dimensionedScalar rhoAir(dict.lookup("rhoAir"));
        
            //rho = rho12 + rhoAir;
            fluxOld = fvc::interpolate(rt.oldTime())*phi;
            flux = fvc::interpolate(rt)*phi;
            
            //S1: Is there evaporation?
            //S2: Is there condensation?
            //S3: Is the evaporation limited by the liquid?
            //S4: Is there enough liquid to reach equilibrium?
            forAll(S1, celli)
            {
                if (rv[celli] >= rvs[celli])
                {
                    S1[celli] = 0;
                    S2[celli] = 1;
                    S3[celli] = 0;
                    S4[celli] = 0;
                }
                else
                {
                    S1[celli] = 1;
                    S2[celli] = 0;
                    
                    if ( (rvs[celli]-rv[celli]) <= rl[celli] )
                    {
                        S3[celli] = 0;
                        S4[celli] = 1;
                    }
                    else
                    {
                        S3[celli] = 1;
                        S4[celli] = 0;
                    }
                }
            }
            Sl = timeScale*S1*S3;
            Sv = timeScale*(S2+S1*S4);
            

            //Prognostic Method.
            rl = rl.oldTime() - dt*
            (
                (1-a)*fvc::div(phi,rl.oldTime())
              +     a*fvc::div(phi,rl_temp)
              + (1-a)*Sl*rl.oldTime()
              - (1-a)*Sv*(rv.oldTime()-rvs)
              -     a*Sv*(rv_temp-rvs)
            );
            if (implicitSource) rl = rl/(dt*a*Sl+1);
            else                rl -= dt*a*Sl*rl_temp;
            
            rv = rv.oldTime() - dt* 
            (
                (1-a)*fvc::div(phi,rv.oldTime())
              +     a*fvc::div(phi,rv_temp)
              - (1-a)*Sl*rl.oldTime()
              + (1-a)*Sv*(rv.oldTime()-rvs)
              -     a*Sv*rvs
            );
            if (implicitSource)
            {
                rv += dt*a*Sl*rl;
                rv =  rv/(dt*a*Sv+1);
                rl -= dt*a*Sv*(rv_temp-rv);
            }
            else                
            {
                rv += dt*a*(Sl*rl_temp - Sv*rv_temp);
            }
            
            rl_temp = rl;
            rv_temp = rv;
            
            
            //Diagnostic Method
            fvScalarMatrix rtEqn
            (
                fvm::ddt(rt)
                + (1-a)*fvc::div(phi, rt.oldTime())
                +     a*fvc::div(phi, rt)
            );
            rtEqn.solve();
            
            rv_diag = min((rt-rhoAir)/rhoAir,rvs);
            rl_diag = (rt-rhoAir)/rhoAir - rv_diag;
        
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
        Info << "Transfer Term Min: " 
             << min(S.internalField()).value() << " Transfer Term Max: "
             << max(S.internalField()).value() << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
