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
    ImplicitAdvectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "velocityField.H"
#include "CourantNoFunc.H"
#include "localMax.H"
#include "fvModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    const int nCorr(readLabel(mesh.solution().lookup("nOuterCorrections")));

    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "velocityFieldDict",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    bool timeVaryingWind = false;
    autoPtr<velocityField> v;
    if(dict.size() > 0)
    {
        Info << "Setting wind from dictionary" << endl;
        timeVaryingWind = bool(dict.lookup("timeVaryingWind"));
        v = velocityField::New(dict);
        v->applyTo(phi);
        //phi.oldTime() = phi;
        U = fvc::reconstruct(phi);
        U.write();
        //phi.write();
    }

    // Initialise up diagnostics
    //scalar TV = totalVariation(T);
    Info << runTime.timeName() << " T goes from " 
         << min(T.internalField()).value() << " to "
         << max(T.internalField()).value() << endl; //" TV = " << TV << endl;
    
    localMax<scalar> maxInterp(mesh);

    while (runTime.run())
    {
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        if (timeVaryingWind)
        {
            Info << "Setting wind field (half way time for 2nd order)" << endl;
            runTime.setTime
            (
                runTime.time().value() - 0.5*runTime.deltaTValue(),
                runTime.timeIndex()
            );
            v->applyTo(phi);
            runTime.setTime
            (
                runTime.time().value() + 0.5*runTime.deltaTValue(),
                runTime.timeIndex()
            );
        }
        #include "CourantNo.H"

//        Co = CourantNo(phi, runTime.deltaT());
//        surfaceScalarField Cof(maxInterp.interpolate(Co));
//        offCentre = max(0.5, 1 - 1/Cof);
        for (int corr = 0; corr < nCorr; corr++)
        {
        
            // Implicit advection
            Info << "Outer corr " << corr << " Implicit/explicit advection"
                 << endl;
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
            );
            Info << "Solving advection equation" << endl;
            TEqn.solve();
        }
        
        //TV = totalVariation(T);
        Info << runTime.timeName() << " T goes from " 
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;//" TV = " << TV << endl;

        if (runTime.writeTime())
        {
            Co = CourantNo(phi, runTime.deltaT());
            runTime.write();
        }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
