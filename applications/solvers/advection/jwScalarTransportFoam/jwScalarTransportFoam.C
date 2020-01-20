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
    Foam::argList::addBoolOption("RK2", "Run RK2 explicitly");
    Foam::argList::addBoolOption("CN", "Run CN implicitly");
    Foam::argList::addBoolOption("FE", "Run FE explicitly");
    Foam::argList::addBoolOption("BE", "Run BE implicitly");
    Foam::argList::addBoolOption("BDF2", "Run BDF2 implicitly");
    #include "setRootCase.H"
    /* J */
    const Switch RK2 = args.options().found("RK2");
    const Switch CN = args.options().found("CN");
    const Switch FE = args.options().found("FE");
    const Switch BE = args.options().found("BE");
    const Switch BDF2 = args.options().found("BDF2");
    
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "jwCreateFields.H"
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = readLabel(itsDict.lookup("nCorr"));
    //Read a scalar stored in the schemes dictionary.
    //const dictionary& schemesDict = mesh.schemesDict();
    //const scalar offCenter = readScalar(schemesDict.lookup("offCenter"));
   // const dictionary& controlDict = mesh.schemesDict();
   // const scalar deltaT = readScalar(schemesDict.lookup("offCenter"));

//   const dictionary& SchemesDict = mesh.schemesDict().subDict("divSchemes");
//    const string lowOrder = readLabel(SchemesDict.lookup("lowOrder"));
  //  const string hiOrder = readLabel(SchemesDict.lookup("hiOrder"));
//    const string correction = readLabel(SchemesDict.lookup("correction")); this isnt needed
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        // Fixed number of iterations per time-step version

        if (CN)
        {
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 0.5*divPhiT.oldTime()/* why not fvc::div(phi,T)*/
	          + 0.5*fvm::div(phi, T)
	          //+ 0.5*fvc::div(phi, T, "highOrder")
	         // - 0.5*fvc::div(phi, T, "lowOrder")*offCenter*1/offCenter
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

        if (FE)
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
        if (BE)
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
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + 1*fvm::div(phi,T)
                );
	        TEqn.solve();
            for (int corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  - 1/(3*runTime.deltaTValue())*T.oldTime()
		  + 1/(3*runTime.deltaTValue())*T.oldTime().oldTime()
		  + (2/3)*divPhiT.oldTime()
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
        //divPhiThi = fvc::div(phi, T, "highOrder");
	divPhiT = fvc::div(phi, T);
	//divPhiT = fvc::div(phi, T, "correction");

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
