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
    implicitAdvectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "linear.H"
#include "surfaceInterpolate.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvmLaplacian.H"
#include "fvmDiv.H"
#include "fvScalarMatrix.H"
#include "velocityField.H"
#include "CourantNoFunc.H"
#include "butcherTableau.H"
#include "RKfield.H"
#include "EulerDdtScheme.H"
#include "localMax.H"
#include "fvModels.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "createFvModels.H"
    #include "createRKFields.H"
    #define aii Bt.coeffs()[iRK][iRK]
    // Pre-defined time stepping scheme and min/max interpolation
    fv::EulerDdtScheme<scalar> EulerDdt(mesh);
    localMax<scalar> maxInterp(mesh);
    const scalar CoLimit = readScalar(mesh.schemes().lookup("CoLimit"));

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
        U = fvc::reconstruct(phi);
        U.write();
    }

    Co = CourantNo(phi, runTime.deltaT());
    Info << "Courant Number mean: "
         << (fvc::domainIntegrate(Co)/Vtot).value()
         << " max: " << max(Co).value() << endl;

    // Implicit/explicit blend
    ImEx = max(scalar(0), 1-CoLimit/maxInterp.interpolate(Co));
    phii = ImEx*phi;
    phie = (1 - ImEx)*phi;
    
    Info << runTime.timeName() << " T goes from " 
         << min(T.internalField()).value() << " to "
         << max(T.internalField()).value() << endl;
    
    while (runTime.run())
    {
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        for (int iRK = 0; iRK < Bt.nSteps(); iRK++)
        {
            if (timeVaryingWind)
            {
                Info << "Setting wind field" << endl;
                runTime.setTime
                (
                    runTime.time().value() - (1-Bt.subTimes()[iRK])*dt,
                    runTime.timeIndex()
                );
                v->applyTo(phi);
                runTime.setTime
                (
                    runTime.time().value() + (1-Bt.subTimes()[iRK])*dt,
                    runTime.timeIndex()
                );
                Co = CourantNo(phi, runTime.deltaT());
                Info << "Courant Number mean: "
                     << (fvc::domainIntegrate(Co)/Vtot).value()
                     << " max: " << max(Co).value() << endl;
            
                if (iRK == 0)
                {
                    ImEx = max(scalar(0), 1-CoLimit/maxInterp.interpolate(Co));
                    phii = ImEx*phi;
                    phie = (1 - ImEx)*phi;
                }
            }

            // Implicit RK stage
            fvScalarMatrix TEqn
            (
                EulerDdt.fvmDdt(T)
              + Bt.subTimes()[iRK]*fvm::div(phii, T)
              + aii*fvm::div(phie, T)
              == Bt.RKstep(iRK, dTdt)
              + fvModels.source(T)
            );
            TEqn.solve();

            // update time derivatives for next stage
            dTdt[iRK] = -fvc::div(phie, T);
        }
        
        if (max(mag(Bt.weights())) > SMALL)
        {
            // Final RK steps
            T = T.oldTime() + dt*Bt.RKfinal(dTdt) - dt*fvc::div(phii, T);
        }
        
        const dimensionedScalar Ttot = fvc::domainIntegrate(T);
        Info << runTime.timeName() << " T goes from " 
             << min(T.internalField()).value() << " to "
             << max(T.internalField()).value()
             << " normalised T mass change = "
             << ((Ttot - Ttot0)/Ttot0).value()
             << endl;

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
