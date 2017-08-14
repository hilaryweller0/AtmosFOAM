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
        
        q1_temp = q1;
        q2_temp = q2;
        
        for (int corr=0; corr < 3; corr++)
        {
            //const dimensionedScalar rhoAir(dict.lookup("rhoAir"));
        
            //rho = rho12 + rhoAir;
            fluxOld = fvc::interpolate(rho.oldTime())*phi;
            flux = fvc::interpolate(rho)*phi;
            
            
            if (BryanFritschSource)
            {
                //S1: Is there evaporation?
                //S2: Is there condensation?
                //S3: Is the evaporation limited by the liquid?
                //S4: Is there enough liquid to reach equilibrium?
                forAll(S1, celli)
                {
                    if (q2[celli] >= rvs[celli])
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
                        
                        if ( (rvs[celli]-q2[celli]) <= q1[celli] )
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
            }
            else
            {
                S = (q2-JahnScale*es/(Rv*T)) - q1 + sqr((q2-JahnScale*es/(Rv*T))*(q2-JahnScale*es/(Rv*T)) + q1*q1);
            }

            // Set up the matrix without adding implicit/explicit parts
            // of advection or source terms
            q1 = q1.oldTime() - dt/(0.5*transferTerm*S1*S3+1) *
            (
                0.5*fvc::div(phi,q1.oldTime())
              + 0.5*fvc::div(phi,q1_temp)
              + 0.5*timeScale*transferTerm*S1*S3*q1.oldTime()
              - 0.5*timeScale*transferTerm*(S2+S1*S4)*(q2.oldTime()-rvs)
              - 0.5*timeScale*transferTerm*(S2+S1*S4)*(q2_temp-rvs)
            );
            
            q2 = q2.oldTime() - dt/(0.5*transferTerm*(S2+S1*S4)+1) * 
            (
                0.5*fvc::div(phi,q2.oldTime())
              + 0.5*fvc::div(phi,q2_temp)
              - 0.5*timeScale*transferTerm*S1*S3*q1.oldTime()
              - 0.5*timeScale*transferTerm*S1*S3*q1
              + 0.5*timeScale*transferTerm*(S2+S1*S4)*(q2.oldTime()-2*rvs)
            );
            
            q1 -= 0.5*dt*timeScale*transferTerm*(S2+S1*S4)*(q2_temp-q2);
            q1_temp = q1;
            q2_temp = q2;
            
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho)
                + fvc::div(phi, rho)
            );

            rhoEqn.solve();
            
            q2_analytic = min((rho-rhoAir)/rhoAir,rvs);
            q1_analytic = (rho-rhoAir)/rhoAir - q2_analytic;
            
            Info << "Writing rho and q after rho.solve" << endl;
            
            
        }
        Info << " rho goes from " << min(rho.internalField()).value() << " to "
             << max(rho.internalField()).value() << endl;
        Info << " q1 goes from " << min(q1.internalField()).value() << " to "
             << max(q1.internalField()).value() << endl;
        Info << " Total rho in system: " << sum(rho.internalField()-1) << endl;
        Info << " Total Moisture: " << sum(q1.internalField())+sum(q2.internalField()) << endl;
        Info << "Transfer Term Min: " << min(S.internalField()).value() << " Transfer Term Max: " << max(S.internalField()).value() << endl;
        runTime.write();
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
