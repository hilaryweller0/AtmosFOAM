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
    ImExAdvectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "velocityField.H"
#include "CourantNoFunc.H"
#include "EulerDdtScheme.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    const int nCorr(readLabel(mesh.solution().lookup("nOuterCorrections")));
    const scalar CoLimit = readScalar(mesh.schemes().lookup("CoLimit"));
    const Switch fullImplicit(mesh.schemes().lookup("fullImplicit"));

    // Pre-defined time stepping scheme and min/max interpolation
    fv::EulerDdtScheme<scalar> EulerDdt(mesh);
    localMax<scalar> maxInterp(mesh);

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
        phi.oldTime() = phi;
        //phi.write();
    }

    #include "CourantNo.H"
    scalar maxCo = CoNum;

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
            Info << "Setting wind field" << endl; // half way time for 2nd order
            runTime.setTime
            (
                runTime.time().value() - 0.5*runTime.deltaTValue(),
                runTime.timeIndex()
            );
            v->applyTo(phi);
            // For variable density flow
            //phiv = phi/rhof;

            #include "CourantNo.H"
            maxCo = CoNum;
            Co = CourantNo(phi, runTime.deltaT());
            Cof = maxInterp.interpolate(Co);
            if (CoLimit > SMALL)
            {
                ImEx = 0.5*(sign(Cof - CoLimit) + 1);
                offCentre = max(0.5, 1 - 1/Cof);
            }

            runTime.setTime
            (
                runTime.time().value() + 0.5*runTime.deltaTValue(),
                runTime.timeIndex()
            );
        }

        Info << "Advection" << endl;
        T -= runTime.deltaT()*fvc::div((1-offCentre)*phi, T, "div(phi,T)");
        T.oldTime() = T;

        for (int corr = 0; corr < nCorr; corr++)
        {
            fvScalarMatrix TEqn
            (
                //fvm::ddt(rho, T)
                EulerDdt.fvmDdt(T)
            );
            if (CoLimit > SMALL)
            {
                TEqn += fvc::div(offCentre*(1-ImEx)*phi, T, "div(phi,T)");
            }
            if (maxCo > CoLimit)
            {
                TEqn += fvScalarMatrix(fvm::div(offCentre*ImEx*phi, T, "div(phi,T)"));
            }
            if (fullImplicit)
            {
                TEqn.solve();
            }
            else
            {
                T = TEqn.H()/TEqn.A();
                //T = T.oldTime()
                //  - runTime.deltaT()*fvc::div(offCentre*(1-ImEx)*phi, T, "div(phi,T)");
            }
        }
        
        //TV = totalVariation(T);
        Info << runTime.timeName() << " T goes from " 
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;//" TV = " << TV << endl;

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
