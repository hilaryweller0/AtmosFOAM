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
#include "upwind.H"
#include "cubicUpwind.H"
#include "gaussConvectionScheme.H"
#include "fvcFluxLimit.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    // Courant number limit for pure explicit/HO
    const scalar CoLimit = readScalar
    (
        mesh.schemes().subDict("divSchemes").lookup("CoLimit")
    );
    // list of the names of the tracers
    const wordList Tnames
    (
        mesh.solution().lookupOrDefault<wordList>("tracerNames", wordList(1, "T"))
    );
    // Is tracer zero density for the other tracers
    const Switch withDensity
    (
        mesh.solution().lookupOrDefault<Switch>("withDensity", false)
    );
    const Switch FCT(mesh.schemes().subDict("divSchemes").lookup("FCT"));
    
    // Pre-defined time stepping scheme and min/max interpolation
    fv::EulerDdtScheme<scalar> EulerDdt(mesh);
    localMax<scalar> maxInterp(mesh);

    Info << "Butcher Tableau\n" << endl;
    const butcherTableau Bt
    (
        mesh.schemes().subDict("ddtSchemes").lookup
        (
            "butcherTableau"
        )
    );

    #include "createFields.H"
    #include "createFvModels.H"
    
    // Low order interpolation and HO correction scheme
    upwind<scalar> lo(mesh, phi);
    fv::gaussConvectionScheme<scalar> upwindConvect
    (
        mesh, phi,
        tmp<surfaceInterpolationScheme<scalar>>(new upwind<scalar>(mesh, phi))
    );
    tmp<surfaceInterpolationScheme<scalar>> hct
    (
        surfaceInterpolationScheme<scalar>::New
        (
            mesh, phi,
            mesh.schemes().subDict("divSchemes").lookup("correctionScheme")
        )
    );
    surfaceInterpolationScheme<scalar>& hc = hct.ref();

    const dimensionedScalar Ttot0 = fvc::domainIntegrate(Tsum);

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary velocityDict
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
    if(velocityDict.size() > 0)
    {
        Info << "Setting wind from dictionary" << endl;
        timeVaryingWind = readBool(velocityDict.lookup("timeVaryingWind"));
        v = velocityField::New(velocityDict);
        v->applyTo(volFlux);
        U = fvc::reconstruct(volFlux);
        U.write();
    }

    Co = CourantNo(volFlux, dt);
    Info << "Courant Number mean: "
         << (fvc::domainIntegrate(Co)/Vtot).value()
         << " max: " << max(Co).value() << endl;

    // Print Diagnostics
    for(label iT = 0; iT < T.size(); iT++)
    {
        Info << T[iT].name() << " goes from "
             << min(T[iT].internalField()).value() << " to "
             << max(T[iT].internalField()).value() << endl;
    }
    Info << runTime.timeName() << " Tsum goes from " 
         << min(Tsum.internalField()).value() << " to "
         << max(Tsum.internalField()).value() << endl;
    
    while (runTime.run())
    {
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        // Loop over each RK stage
        for (int iRK = 0; iRK < Bt.nSteps(); iRK++)
        {
            const scalar ci = Bt.subTimes()[iRK];
            Info << "iRK = " << iRK << endl;
            if (timeVaryingWind)
            {
                Info << "Setting wind field" << endl;
                runTime.setTime
                (
                    runTime.time().value() - (1-ci)*dt.value(),
                    runTime.timeIndex()
                );
                v->applyTo(volFlux);
                runTime.setTime
                (
                    runTime.time().value() + (1-ci)*dt.value(),
                    runTime.timeIndex()
                );
                Co = CourantNo(volFlux, dt);
                Info << "Courant Number mean: "
                     << (fvc::domainIntegrate(Co)/Vtot).value()
                     << " max: " << max(Co).value() << endl;
                if (iRK == 0)
                {
                    beta = max(scalar(0), 1-CoLimit/maxInterp.interpolate(Co));
                    alpha = max
                    (
                        scalar(0.5), //beta
                        1-2/maxInterp.interpolate(Co)
                    );
                }
            }
            // Advection for each tracer
            for(int iT = 0; iT < T.size(); iT++)
            {
                if (iRK == 0)
                {
                    fluxOld[iT] = (1-alpha)*beta*phi*lo.interpolate(T[iT]);
                    if (FCT && !(withDensity && iT == 0))
                    {
                        fluxOld[iT].oldTime() = (1-alpha*beta)*phi
                                                *lo.interpolate(T[iT]);
                    }
                }
                
                if (mag(ci) > SMALL)
                {
                    fluxT[iT] = ci*fluxOld[iT];
                    for(int lRK = 0; lRK < iRK; lRK++)
                    {
                        fluxT[iT] += Bt.coeffs()[iRK][lRK]*fluxTmp[iT][lRK];
                    }

                    // Implicit RK stage
                    fvScalarMatrix TEqn
                    (
                        upwindConvect.fvmDiv(alpha*beta*ci*phi, T[iT])
                    );
                    if (withDensity && iT >= 1)
                    {
                        TEqn += fvScalarMatrix(EulerDdt.fvmDdt(T[0], T[iT]))
                             - ci*fvModels.source(T[0], T[iT]);
                    }
                    else
                    {
                        TEqn += fvScalarMatrix(EulerDdt.fvmDdt(T[iT]))
                             - ci*fvModels.source(T[iT]);
                    }
                    
                    // Low order solution (for FCT)
                    if (FCT && !(withDensity && iT == 0))
                    {
                        fluxT[iT].oldTime() = ci*fluxOld[iT].oldTime();
                        solve
                        (
                            TEqn == -fvc::div(fluxT[iT].oldTime())
                        );
                        Info << "Low order " << T[iT].name() << " goes from "
                             << min(T[iT].field()) << " to "
                             << max(T[iT].field()) << endl;
                        // Store low order solution and low order flux
                        T[iT].oldTime().oldTime() = T[iT];
                        fluxT[iT].oldTime() += TEqn.flux();
                    }
                    
                    // HO solution
                    solve(TEqn == -fvc::div(fluxT[iT]));
                    Info << "Full " << T[iT].name() << " goes from "
                         << min(T[iT].field()) << " to "<< max(T[iT].field())
                         << endl;
                    
                    fluxT[iT] += TEqn.flux();
                } // end if (mag(ci) > SMALL)
                
                // Full HO fluxes
                fluxTmp[iT][iRK] = phi*(1-beta)*lo.interpolate(T[iT]);
                if (withDensity && iT >= 1)
                {
                    fluxTmp[iT][iRK] += volFlux*
                    (
                        hc.interpolate(T[0])*hc.correction(T[iT])
                      + hc.correction(T[0])*lo.interpolate(T[iT])
                    );
                }
                else {fluxTmp[iT][iRK] += phi*hc.correction(T[iT]);}
                
                // Modify the fluxes and solution for FCT
                if (FCT && !(withDensity && iT == 0))
                {
                    surfaceScalarField fluxCorr = fluxT[iT] - fluxT[iT].oldTime();

                    fvc::fluxLimit(fluxCorr, T[0]*T[iT].oldTime().oldTime(), dt);

                    fluxT[iT] = fluxT[iT].oldTime() + fluxCorr;
                    if (withDensity && iT > 0)
                    {
                        //T[iT] =  (T[0].oldTime()*T[iT].oldTime() 
                        //         - dt*fvc::div(fluxT[iT]))/T[0];
                        T[iT] = T[iT].oldTime().oldTime() 
                              - dt*fvc::div(fluxCorr)/T[0];
                    }
                    else
                    {
                        T[iT] =  T[iT].oldTime() - dt*fvc::div(fluxT[iT]);
                    }
                    Info << "Limited " << T[iT].name() << " goes from "
                         << min(T[iT].field()) << " to "
                         << max(T[iT].field()) << endl;
                } // end of FCT
                
                // Replace volume flux with mass flux if withDensity
                if (withDensity && iT == 0)
                {
                    massFlux = phi*lo.interpolate(T[0]);
                    phip = &massFlux;
                }
            } // End of Advection for each tracer
            if (withDensity)
            {
                phip = &volFlux;
            }
        } // End of Loop over each RK stage
        
        // Print Diagnostics
        for(label iT = 0; iT < T.size(); iT++)
        {
            Info << T[iT].name() << " goes from " << min(T[iT].field())
                 << " to " << max(T[iT].field()) << endl;
        }
        
        // Sum of all tracers
        if (T.size() > 1)
        {
            if (withDensity) {Tsum = T[0]*T[1];}
            else
            {
                Tsum = T[0];
                for(label iT = 1; iT < T.size(); iT++)
                {
                    Tsum += T[iT];
                }
            }
        }
        const volScalarField& Ttmp = T.size() > 1 ? Tsum : T[0];
        const dimensionedScalar Ttot = fvc::domainIntegrate(Ttmp);
        Info << runTime.timeName() << " Tsum goes from " 
             << min(Ttmp.internalField()).value() << " to "
             << max(Ttmp.internalField()).value()
             << " normalised Tsum mass change = "
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
