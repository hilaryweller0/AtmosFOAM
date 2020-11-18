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
    scalarDeformation

Description
    Solves a transport equation for a passive scalar using Implicit
    time-stepping for the evolving velocity field on a plane
    for deformational flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "deformationalFlow.H"
#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "fvcFluxLimit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = readLabel(itsDict.lookup("nCorr"));
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    const Switch implicit = mesh.schemesDict().lookup("implicit");
    const dictionary& divSchemesDict = mesh.schemesDict().subDict("divSchemes");
    const Switch FCT = mesh.schemesDict().lookup("fluxCorrectedTransport");

    // Create the class for the deformational flow
    deformationalFlow defFlow
    (
        IOdictionary
        (
            IOobject
            (
                "deformationalAdvectionDict", "system", runTime,
                IOobject::MUST_READ
            )
        )
    );

    #include "createFields.H"

    // Pre-define space schemes
    upwind<scalar> upwindInterp(mesh, phi);
    fv::gaussConvectionScheme<scalar> upwind(mesh, phi, upwindInterp);
    fv::gaussConvectionScheme<scalar> divScheme
    (
        mesh, phi, divSchemesDict.lookup("div(phi,T)Interpolation")
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        defFlow.update(phi, U, Uf);
        #include "CourantNo.H"

        // Mid-time T calculating
        surfaceScalarField timeFlux = dt*(1-offCentre)*phi.oldTime();
        // High order flux
        surfaceScalarField fluxCorr
             = timeFlux*divScheme.interpolate(timeFlux, T.oldTime());
        if (FCT && implicit)
        {
            // Low order update
            fvScalarMatrix TEqn
            (
                fvm::Sp(1, T) - T.oldTime()
              + upwind.fvmDiv(timeFlux, T)
            );
            TEqn.solve();

            // High order correction
            fluxCorr -= timeFlux*upwind.interpolate(timeFlux, T);
        }
        else if (FCT)
        {
            // Low-order update
            T = T.oldTime() - upwind.fvcDiv(timeFlux, T.oldTime());
            // High-order correction
            fluxCorr -= timeFlux*upwind.interpolate(timeFlux, T.oldTime());
        }
        
        // Limit and apply the correction
        if (FCT) fvc::fluxLimit(fluxCorr, T, T.oldTime());
        T -= fvc::div(fluxCorr);
        
        // New T and timeFlux for the rest of the time step
        T.oldTime() = T;
        timeFlux = dt*offCentre*phi;

        for (int corr = 0; corr < nCorr; corr++)
        {
            // High order flux
            fluxCorr = timeFlux*divScheme.interpolate(timeFlux,T);

            if (implicit)
            {
                // Low order update
                fvScalarMatrix TEqn
                (
                    fvm::Sp(1, T) - T.oldTime()
                  + upwind.fvmDiv(timeFlux, T)
                );
                TEqn.solve();

                // High order correction
                fluxCorr -= timeFlux*upwind.interpolate(timeFlux, T);
            }
            else
            {
                // Low order update
                T = T.oldTime() - upwind.fvcDiv(timeFlux, T.oldTime());

                // High order correction
                fluxCorr -= timeFlux*upwind.interpolate(timeFlux, T.oldTime());
            }
            // Limit and apply the correction
            if (FCT) fvc::fluxLimit(fluxCorr, T, T.oldTime());
            T -= fvc::div(fluxCorr);
        }

        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << nl << endl;
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
