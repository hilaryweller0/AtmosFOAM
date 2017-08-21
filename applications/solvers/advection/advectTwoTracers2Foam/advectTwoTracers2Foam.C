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
        
        rl_temp = rl;
        rv_temp = rv;
        
        for (int corr=0; corr < 3; corr++)
        {
            //const dimensionedScalar rhoAir(dict.lookup("rhoAir"));
        
            //rho = rho12 + rhoAir;
            fluxOld = fvc::interpolate(rt.oldTime())*phi;
            flux = fvc::interpolate(rt)*phi;
            
            
            if (BryanFritschSource)
            {
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
                Sl = timeScale*transferTerm*S1*S3;
                Sv = timeScale*transferTerm*(S2+S1*S4);
            }
            else
            {
                S = (rv-JahnScale*es/(Rv*T)) - rl + sqr((rv-JahnScale*es/(Rv*T))*(rv-JahnScale*es/(Rv*T)) + rl*rl);
            }

            // Set up the matrix without adding implicit/explicit parts
            // of advection or source terms
            rl = rl.oldTime() - dt*
            (
                0.5*fvc::div(phi,rl.oldTime())
              + 0.5*fvc::div(phi,rl_temp)
              + 0.5*Sl*rl.oldTime()
              - 0.5*Sv*(rv.oldTime()-rvs)
              - 0.5*Sv*(rv_temp-rvs)
            );
            rl = rl/(0.5*dt*Sl+1);
            
            rv = rv.oldTime() - dt* 
            (
                0.5*fvc::div(phi,rv.oldTime())
              + 0.5*fvc::div(phi,rv_temp)
              - 0.5*Sl*rl.oldTime()
              - 0.5*Sl*rl
              + 0.5*Sv*(rv.oldTime()-2*rvs)
            );
            rv = rv/(0.5*dt*Sv+1);
            
            rl -= 0.5*dt*Sv*(rv_temp-rv);
            rl_temp = rl;
            rv_temp = rv;
            
            fvScalarMatrix rtEqn
            (
                fvm::ddt(rt)
                + fvc::div(phi, rt)
            );

            rtEqn.solve();
            
            rv_diag = min((rt-rhoAir)/rhoAir,rvs);
            rl_diag = (rt-rhoAir)/rhoAir - rv_diag;
            
            Info << "Writing r rt.solve" << endl;
        }
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
