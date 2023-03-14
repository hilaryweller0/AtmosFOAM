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
    Solves a transport equation for a passive scalar using an explicit
    time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "velocityField.H"
#include "CourandNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("leapfrog", "leapfrog timestepping");
    Foam::argList::addBoolOption("rk2", "2nd order rk time stepping with nCorr iterations. Use nCorr = 2 or 3 (in fvSolution)");
    Foam::argList::addBoolOption("rk3", "Third-order Runge-Kutta scheme used by Skamarock & Gassmann 2011");
    Foam::argList::addBoolOption("rk4", "the classical Runge-Kutta scheme");
    Foam::argList::addBoolOption("forwardEuler", "");
    Foam::argList::addBoolOption("forwardBackward", "");
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    Foam::argList::addBoolOption("explicitTimestepping", "halts if Co > 1");
    Foam::argList::addBoolOption("backwardEuler", "Euler implicit");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initEnergy.H"
    #include "energy.H"

    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solution().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<label>("nCorr", label(2));

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "advectionDict",
            mesh.time().system(),
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

    scalar maxCoNum = 0;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        #include "CourantNo.H"

        if (CoNum > maxCoNum) maxCoNum = CoNum;
        if (explicitTimestepping && CoNum > 1.0)
        {
            FatalErrorInFunction << "Max Courant number > 1"
                 << exit(FatalError);
        }

        if (args.options().found("forwardEuler"))
        {
            T = T.oldTime() - dt*fvc::div(phi,T);
        }
        else if (args.options().found("forwardBackward"))
        {
            T = T.oldTime() - dt*fvc::div(phi,T);
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
        else if (args.options().found("rk4"))
        {
            k1 = -fvc::div(phi, T.oldTime(), "div(phi,T)");
            k2 = -fvc::div(phi, T.oldTime() + dt/2 * k1, "div(phi,T)");
            k3 = -fvc::div(phi, T.oldTime() + dt/2 * k2, "div(phi,T)");
            k4 = -fvc::div(phi, T.oldTime() + dt * k3, "div(phi,T)");
            T = T.oldTime() + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        }
        else if (args.options().found("backwardEuler"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + fvm::div(phi, T)
                );
                TEqn.solve();
            }
        }
        else if (args.options().found("rk2"))
        {
            for (int corr=0; corr < nCorr; corr++)
            {
                T = T.oldTime() - 0.5*dt*
                (
                    fvc::div(phi, T) + fvc::div(phi, T.oldTime())
                );
            }
        }
        else
        {
            FatalErrorIn("advectionFoam") 
                << " no recognised time stepping scheme in options. "
                << "See advectionFoam -help for valid options" 
                << exit(FatalError);
        }

        T.correctBoundaryConditions();
        
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        runTime.write();

        #include "energy.H"

        if (timeVaryingWind)
        {
            v->applyTo(phi);

            U = fvc::reconstruct(phi);
            Uf = linearInterpolate(U);
            Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        }
    }

    Info << "Maximum Courant Number at any timestep: " << maxCoNum << endl;
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
