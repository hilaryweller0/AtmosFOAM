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
#include "velocityField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("leapfrog", "use leapfrog timestepping scheme rather than RK2");
    Foam::argList::addBoolOption("rk2", "two-stage second-order Runge-Kutta");
    Foam::argList::addBoolOption("rk3", "Third-order Runge-Kutta scheme used by Skamarock & Gassmann 2011");
    Foam::argList::addBoolOption("rk4", "the classical Runge-Kutta scheme");
    Foam::argList::addBoolOption("forwardEuler", "");
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    Foam::argList::addBoolOption("explicitTimestepping", "halts if Co > 1");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "advectionDict",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    bool explicitTimestepping = args.options().found("explicitTimestepping");
    bool timeVaryingWind = dict.lookupOrDefault<bool>("timeVaryingWind", false);
    const dictionary& velocityDict = dict.subOrEmptyDict("velocity");
    autoPtr<velocityField> v;
    if (velocityDict.size() > 0)
    {
        v = velocityField::New(dict.subOrEmptyDict("velocity"));
    }
    
    // go backwards in time by one time step to initialise leap-frog
    T.oldTime().oldTime() = T + dt*fvc::div(phi, T);

    while (runTime.loop())
    {
        #include "CourantNo.H"
        if (explicitTimestepping && CoNum > 1.0)
        {
            FatalErrorInFunction << "Max Courant number > 1" << exit(FatalError);
        }
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (args.options().found("forwardEuler"))
        {
            T = T.oldTime() - dt*fvc::div(phi,T);
        }
        else if (args.options().found("leapfrog"))
        {
            T = T.oldTime().oldTime() - 2*dt*fvc::div(phi,T);
        }
        else if (args.options().found("rk3"))
        {
            T = T.oldTime() - dt/3*fvc::div(phi, T);
            T = T.oldTime() - dt/2*fvc::div(phi, T);
            T = T.oldTime() - dt*fvc::div(phi, T);
        }
        else if (args.options().found("rk2"))
        {
            T = T.oldTime() - 0.5*dt*fvc::div(phi,T.oldTime());
            T = T.oldTime() - dt*fvc::div(phi,T);
        }
        else if (args.options().found("rk4"))
        {
            k1 = -fvc::div(phi, T.oldTime(), "div(phi,T)");
            k2 = -fvc::div(phi, T.oldTime() + dt/2 * k1, "div(phi,T)");
            k3 = -fvc::div(phi, T.oldTime() + dt/2 * k2, "div(phi,T)");
            k4 = -fvc::div(phi, T.oldTime() + dt * k3, "div(phi,T)");
            T = T.oldTime() + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
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
            v->applyTo(phi);

            U = fvc::reconstruct(phi);
            Uf = linearInterpolate(U);
            Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
