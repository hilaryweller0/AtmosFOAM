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
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar using Implicit
    time-stepping. Option -explicit runs everything explicitly

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CourantNo.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = readLabel(itsDict.lookup("nCorr"));
    
    // The off-centering and the choice of anti-diffusive flux
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    CourantNo(phi, dt);

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;

        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + (1-offCentre)*fvc::div(phi.oldTime(), T, "upwind")
          + offCentre*fvm::div(phi, T, "upwind")
        );
        TEqn.solve();

        T.oldTime() = T;

        // Fixed number of iterations per time-step version
        for (int corr = 0; corr < nCorr; corr++)
        {
            surfaceVectorField gradT = linearInterpolate(fvc::grad(T));
            //gradT += (fvc::snGrad(T) - (gradT & mesh.Sf())/mesh.magSf())
            //        *mesh.Sf()/mesh.magSf();

            surfaceScalarField Tf = stabilise
            (
                linearInterpolate(T),
                dimensionedScalar("", T.dimensions(), SMALL)
            );

            antiD = 0.5/Tf*
            (
                mag(phi)*fvc::snGrad(T)/mesh.deltaCoeffs()
              - (1-2*offCentre)*phi*dt* (Uf & gradT)
            );
            
            CourantNo(antiD, dt);

            // Limit the anti-diffusive velocity so that Courant<1
            antiD = min(antiD, 0.5*mesh.magSf()/mesh.deltaCoeffs()/dt);
            antiD = max(antiD, -0.5*mesh.magSf()/mesh.deltaCoeffs()/dt);

            CourantNo(antiD, dt);
    
            // Ante-diffusive velocity to write out
            antiDv = linearInterpolate(fvc::reconstruct(antiD));
            antiDv += (antiD - (antiDv & mesh.Sf()))*mesh.Sf()
                /sqr(mesh.magSf());
            
            // Apply the correction explicitly
            T = T.oldTime() - dt*fvc::div(antiD, T, "upwind");

            Info << "T goes from " << min(T).value() << " to "
                 << max(T).value() << endl;
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
