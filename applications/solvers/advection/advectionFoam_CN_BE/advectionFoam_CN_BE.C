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
    advectionFoam_CN_BE

Description
    Solves a transport equation for a passive scalar using Crank-Nicolson
    which blends to Backard Euler for large Courant number

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CourantNoFunc.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    // Calculate Courant number
    volScalarField CourantC = CourantNo(phi, dt);
    CourantC.rename("CourantC");
    CourantC.write();
    Info << "Courant number goes from " << min(CourantC).value() << " to "
         << max(CourantC).value() << endl;
    
    // Calculate off centering on the face
    volScalarField offCentreC
    (
        IOobject("offCentre", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, 0.5)
    );
    offCentreC = max(offCentreC, 1 - 1/(CourantC + SMALL));
    offCentreC.write();
    Info << "Off centering goes from " << min(offCentreC).value() << " to "
         << max(offCentreC).value() << endl;
    surfaceScalarField offCentre = linearInterpolate(offCentreC);
    
    // Velocity field divergence free?
    volScalarField divu
    (
        "divu",
        fvc::div(offCentre*phi + (1-offCentre)*phi.oldTime())
    );
    Info << "divu goes from " << min(divu).value() << " to "
         << max(divu).value() << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvc::div((1-offCentre)*phi.oldTime(), T.oldTime(), "upwind")
          + fvm::div(offCentre*phi, T, "upwind")
        );
        TEqn.solve();

        Info << "T goes from " << min(T).value() << " to "
             << max(T).value() << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
