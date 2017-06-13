/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    advectionConservativeF

Description
    Conservative advection on faces using fvc::interpolate(fvc::div())

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volInterpolationScheme.H"
#include "gravity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createGravity.H"
    #include "createVolInterpolation.H"
    #include "createFields.H"

    Info << "\nCalculating advection\n" << endl;

    #include "CourantNo.H"

    Tf = Tf * mag(g.unitFaceNormal()) + (1.0 - mag(g.unitFaceNormal()))*fvc::interpolate(T);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int corr=0; corr < 3; corr++)
        {
            Tf = Tf.oldTime() - 0.5*dt *
            (
                (fvc::interpolate(fvc::div(Tf*phi))) + 
                (fvc::interpolate(fvc::div(Tf.oldTime()*phi.oldTime())))
            );
        }

        T = cellCentreReconstruction.interpolate(Tf);

        Tf = Tf * mag(g.unitFaceNormal()) + (1.0 - mag(g.unitFaceNormal()))*fvc::interpolate(T);

        Info << " Tf goes from " << min(Tf.internalField()) << " to "
             << max(Tf.internalField()) << endl;
        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << endl;
        runTime.write();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
