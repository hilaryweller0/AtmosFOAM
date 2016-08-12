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
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;


    bool firstStep = true;

    while (runTime.run())
    {
        const volVectorField U
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

        const surfaceVectorField Uf
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

        const surfaceScalarField phi
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

        if (firstStep)
        {
            // go backwards in time by one time step to initialise leap-frog
            T.oldTime().oldTime() = T + dt*fvc::div(phi, T);
            firstStep = false;
        }

        #include "CourantNo.H"

        runTime.loop();
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
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
