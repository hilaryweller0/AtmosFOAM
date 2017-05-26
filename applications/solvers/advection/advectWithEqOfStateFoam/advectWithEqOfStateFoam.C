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
    advectWithEqOfStateFoam

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
            fvScalarMatrix rhoEqn(fvm::ddt(rho) + 0.5*fvc::div(phi, rho.oldTime()));

            // Add the advection terms either implicit or explicit
            if (implicitAdvection)
            {
                rhoEqn += 0.5*fvm::div(phi, rho);
            }
            else
            {
                rhoEqn += 0.5*fvc::div(phi, rho);
            }
            
            
            // Solve the matrices for the equations
            rhoEqn.solve();
        }
        
        double t_now = runTime.time().value();
        double angVel = -0.005236;
        Info << 1 << endl;
        //volScalarField q = 0*r/(r+sourceLengthScale);
        Info << 2 << endl;
        
        if (t_now <= 300.)
        {
            //q = Foam::exp(r/sourceLengthScale * (Foam::sin(angVel*t_now + 0.5*M_PI)-1) );
        }        
        else 
        {
            //q = 1 - (1 - Foam::exp(-2*r/sourceLengthScale))*Foam::exp(-r/sourceLengthScale * (Foam::sin(angVel*t_now + 0.5*M_PI)+1) );
        }
        
        //for (dimensionedScalar loop=0*runTime.timeName()*w; loop+1 < runTime.timeName()*w/(2*M_PI); loop = loop + 1)
        //{}
        
        rho1 = rho;
        
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
