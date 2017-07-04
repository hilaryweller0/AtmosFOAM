/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    partitionedShallowWaterFoam

Description
    Transient Solver for shallow water partitioned flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 2);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int icorr=0; icorr < nCorr; icorr++)
        {
            // Solve the continuity equation for each partition and sum
            h == dimensionedScalar("h", dimLength, scalar(0));
            for(label ip = 0; ip < nParts; ip++)
            {
                hSigma[ip] = hSigma[ip] - dt*fvc::div(flux[ip]);
                h += hSigma[ip];
            }
            // Calculate sigma as a diagnostic
            for(label ip = 0; ip < nParts; ip++)
            {
                sigma[ip] = hSigma[ip]/h;
                Info << "Minimim " << sigma[ip].name() << " = "
                     << min(sigma[ip]).value() << endl;
            }

            // Solve the momentum equation for each partition to update the
            // flux and u in each partition
            for(label ip = 0; ip < nParts; ip++)
            {
                flux[ip] = flux[ip].oldTime() - dt*
                (
                    (fvc::interpolate(fvc::div(flux[ip], u[ip])) & mesh.Sf())
                  + g*fvc::interpolate(hSigma[ip])*fvc::snGrad(h)*mesh.magSf()
                );
                u[ip] = fvc::reconstruct(flux[ip])/hSigma[ip];
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
