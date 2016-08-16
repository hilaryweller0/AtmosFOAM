/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    moistThermoVars

Description
    Calculates moist thermodynamic variables from Exner and thetaRho
    and from moisture variables qv and ql

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "atmosphere.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
        atmos.readUpdate();

        Info << "Reading theta" << endl;
        volScalarField theta
        (
            IOobject("theta", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        volScalarField Exner("Exner", atmos.ExnerFromTheta(theta));
        Exner.write();
        volScalarField T("T", Exner*theta);
        T.write();
        volScalarField p("p", air.pFromExner(Exner));
        p.write();
        volScalarField thetae("thetae", atmos.thetae(T));
        thetae.write();
        volScalarField rho("rho", atmos.sumDensity());
        rho.write();
        volScalarField rv("rv", water.gas().rho()/air.rho());
        rv.write();
        volScalarField rl("rl", water.liquid().v()*water.liquid().rho()
                                /air.rho());
        rl.write();
    }

    Info << endl;

    return(0);
}


// ************************************************************************* //
