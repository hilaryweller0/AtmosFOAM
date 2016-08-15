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
    advectionFoam

Description
    Solves a transport equation for a passive scalar using explicit leap-frog
    time-stepping or RK2

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("leapfrog", "use leapfrog timestepping scheme rather than RK2");
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    bool timeVaryingWind = args.options().found("timeVaryingWind");
    
    // go backwards in time by one time step to initialise leap-frog
    T.oldTime().oldTime() = T + dt*fvc::div(phi, T);

    while (runTime.loop())
    {
        #include "CourantNo.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (args.options().found("leapfrog"))
        {
            T = T.oldTime().oldTime() - 2*dt*fvc::div(phi,T);
        }
        else
        {
            for (int corr=0; corr < 3; corr++)
            {
                T = T.oldTime() - 0.5*dt*
                (
                    fvc::div(phi, T) + fvc::div(phi, T.oldTime())
                );
            }
        }

        T.correctBoundaryConditions();
        
        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << endl;
        runTime.write();

        if (timeVaryingWind)
        {
            U = volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("U", dimVelocity, vector::zero),
                "zeroGradient"
            );

            Uf = surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                linearInterpolate(U)
            );

            phi = surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                Uf & mesh.Sf()
            );
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
