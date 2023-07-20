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
    implicitAdvectionFoam1

Description
    Solves a transport equation for a passive scalar using Backward-Euler
    upwind wit H/A instead of a solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "velocityField.H"
#include "CourantNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    const Switch fullSolver(mesh.schemes().lookup("fullSolver"));

    #include "createFields.H"
    const dimensionedScalar Ttot0 = fvc::domainIntegrate(T);

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
        timeVaryingWind = readBool(dict.lookup("timeVaryingWind"));
        v = velocityField::New(dict);
        v->applyTo(phi);
        U = fvc::reconstruct(phi/rhof);
        U.write();
    }

    // Initialise up diagnostics
    //scalar TV = totalVariation(T);
    Info << runTime.timeName() << " T goes from " 
         << min(T.internalField()).value() << " to "
         << max(T.internalField()).value() << endl; //" TV = " << TV << endl;
    
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
            phi.oldTime() = phi;
            phiv = phi/rhof;
            runTime.setTime
            (
                runTime.time().value() + 0.5*runTime.deltaTValue(),
                runTime.timeIndex()
            );
        }
        Co = CourantNo(phiv, runTime.deltaT());
        Info << "Courant Number mean: " << (fvc::domainIntegrate(Co)/Vtot).value()
             << " max: " << max(Co).value() << endl;

        fvScalarMatrix TEqn
        (
            fvm::ddt(rho, T)
          + fvm::div(phi, T)
        );
        if (fullSolver)
        {
            TEqn.solve();
            T = T.oldTime() - runTime.deltaT()*fvc::div(TEqn.flux());
        }
        else
        {
            T = TEqn.H()/TEqn.A();
            T = T.oldTime() - runTime.deltaT()*fvc::div(phi, T);
        }
        
        dimensionedScalar Ttot = fvc::domainIntegrate(T);

        //TV = totalVariation(T);
        Info << runTime.timeName() << " T goes from " 
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value()
             << " normalised T mass change = " << ((Ttot - Ttot0)/Ttot0).value()
             << endl;//" TV = " << TV << endl;

        if (runTime.writeTime())
        {
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
