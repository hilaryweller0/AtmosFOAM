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
    EulerImplicitExplicitAdvectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "velocityField.H"
#include "CourantNoFunc.H"
#include "fvModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvModels.H"
    #include "createTimeControls.H"

    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "velocityFieldDict",
            runTime.system(),
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
        U = fvc::reconstruct(phi);
        U.write();
        phi.write();
    }

    // Initialise up diagnostics
    //scalar TV = totalVariation(T);
    Info << runTime.timeName() << " T goes from " 
         << min(T.internalField()).value() << " to "
         << max(T.internalField()).value() << endl; //" TV = " << TV << endl;
    
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        fvModels.correct();
        #include "CourantNo.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        if (timeVaryingWind)
        {
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

            U = fvc::reconstruct(phi);
//            Uf = linearInterpolate(U);
//            Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        }

        // Implicit advection
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
        );
        TEqn.solve();
        
        //TV = totalVariation(T);
        Info << runTime.timeName() << " T goes from " 
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;//" TV = " << TV << endl;

        if (runTime.writeTime()) {Co = CourantNo(phi, runTime.deltaT());}
        runTime.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
