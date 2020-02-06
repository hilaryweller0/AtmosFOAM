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
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar using Implicit
    time-stepping. Option -explicit runs everything explicitly

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "jwCreateFields.H"
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = readLabel(itsDict.lookup("nCorr"));
    //Read a scalar stored in the schemes dictionary.
    const dictionary& schemesDict = mesh.schemesDict();
    const string timeScheme = schemesDict.lookup("timeScheme");
    const scalar offCenter = readScalar(schemesDict.lookup("offCenter"));
   // const dictionary& controlDict = mesh.schemesDict();
   // const scalar deltaT = readScalar(schemesDict.lookup("offCenter"));

//   const dictionary& SchemesDict = mesh.schemesDict().subDict("divSchemes");
//    const string lowOrder = readLabel(SchemesDict.lookup("lowOrder"));
  //  const string hiOrder = readLabel(SchemesDict.lookup("hiOrder"));
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;
    #include "CourantNo.H"



if (timeScheme == "MPDATA")     // if (timeScheme == "FElowsweeps")
{
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvc::div(phi, T, "lowOrder")
        );
        TEqn.solve();
        
        for (int corr = 0; corr < nCorr; corr++)
        {
            surfaceScalarField antiD = 0.5*fvc::snGrad(T)/linearInterpolate(T)*
            (//the above is incorrect as the derivative should be in the x direction
                mag(phi)/mesh.deltaCoeffs() - sqr(phi)*runTime.deltaT()/mesh.magSf()
            );
            
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvc::div(antiD, T, "lowOrder")
            );
            TEqn.solve();
	    }

        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << nl << endl;

        runTime.write();
    }
}

if (FElowsweeps)     // if (timeScheme == "FElowsweeps")
{
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;
        for (int corr = 0; corr < nCorr; corr++)
        {
	    T = 0;
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + divPhiTlo.oldTime()/10
            );
            TEqn.solve();
            T.oldTime() = T;
	    }
	    divPhiThi = fvc::div(phi, T, "highOrder");
   	    divPhiTlo = fvc::div(phi, T, "lowOrder");
        divPhiT = fvc::div(phi, T);
        divPhiTcorr = fvc::div(phi, T, "correction");

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;

        runTime.write();
    }
}

if (CNlowsweeps)
{
    while (runTime.loop())
    {Info<< "Time = " << runTime.timeName() << endl;

        for (int corr = 0; corr < nCorr; corr++)
        {


            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + 0.5*divPhiTlo.oldTime()
	      + 0.5*fvm::div(phi, T, "lowOrder")


            );
            TEqn.solve();
        T.oldTime() = T;/* this doesnt seem to work */
	}
    divPhiThi = fvc::div(phi, T, "highOrder");
    divPhiTlo = fvc::div(phi, T, "lowOrder");
    divPhiT = fvc::div(phi, T);
    divPhiTcorr = fvc::div(phi, T, "correction");

    Info << " T goes from " << min(T.internalField()) << " to "
         << max(T.internalField()) << nl << endl;
    runTime.write();
    }
}

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        // Fixed number of iterations per time-step version

        if (timeScheme =="CNlow")
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 0.5*divPhiTlo.oldTime()
	          + 0.5*fvm::div(phi, T, "lowOrder")


                );
                TEqn.solve();
            }
        }

        if (CNcorr)
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 0.5*divPhiTcorr.oldTime()
	          + 0.5*fvm::div(phi, T, "correction")
  			

                );
                TEqn.solve();
            }
        }
        if (timeScheme == "CNhigh")
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                ( 
                    fvm::ddt(T)
                  + 0.5*0.5*divPhiTlo.oldTime()
                  + 0.5*0.5*divPhiTcorr.oldTime()
		
	          + 0.5*0.5*fvm::div(phi, T, "lowOrder")
	          + 0.5*0.5*fvm::div(phi, T, "correction")
	
                );
                TEqn.solve();
		cout << "CNhigh implemented";
            }
        }
        if (CN_implicit_low_explicit_corr)
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 0.5*divPhiThi.oldTime()
	              + 0.5*fvm::div(phi, T, "lowOrder")
	              + 0.5*fvc::div(phi, T, "highOrder")
	              - 0.5*fvc::div(phi, T, "lowOrder")
                );
                TEqn.solve();
            }
        }
        if (RK2)
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 0.5*divPhiT.oldTime()
	          + 0.5*fvc::div(phi, T)
                );
	        TEqn.solve();
            }
	}

        if (timeScheme == "FE")
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 1*divPhiT.oldTime()
                );
	        TEqn.solve();
            }
	}
	
        if (timeScheme == "BE")
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 1*fvm::div(phi,T)
                );
	        TEqn.solve();
            }
	}
        if (BDF2)
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  - 1/(3*runTime.deltaT())*T.oldTime()
		  + 1/(3*runTime.deltaT())*T.oldTime().oldTime()
                  + 2/3.*fvm::div(phi,T)
                );
	        TEqn.solve();
            }
	}
        Info << "Max T = " << max(T) << " min T = " << min(T) << endl;

/*        // Keep doing iterations until converged version
        for(bool converged = false; !converged;)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + 0.5*fvm::div(phi, T)
              + 0.5*divPhiT.oldTime()
            );
            solverPerformance sp = TEqn.solve();
            converged = sp.nIterations() <= 1;
        }
*/    
        divPhiThi = fvc::div(phi, T, "highOrder");
        divPhiTlo = fvc::div(phi, T, "lowOrder");
        divPhiT = fvc::div(phi, T);
        divPhiTcorr = fvc::div(phi, T, "correction");

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
