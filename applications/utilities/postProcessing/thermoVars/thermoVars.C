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
    thermoVars

Description
    Calculates thermodynamic variables from Exner and theta

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"

    Info<< "Creating field g.h at cell centres\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        //#   include "createMesh.H"
        mesh.readUpdate();
        gh = g & mesh.C();

        Info << "Reading Exner function" << endl;

        volScalarField Exner
        (
            IOobject("Exner", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        Info << "Reading theta" << endl;

        volScalarField theta
        (
            IOobject("theta", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        Info << "Calculating pressure, p" << flush;
        volScalarField p("p", pRef*pow(Exner, 1./kappa));
        Info << " goes from " << min(p).value() << " to " << max(p).value()
             << endl;
        p.write();

        Info << "Calculating temperature, T" << flush;
        volScalarField T("T", theta*Exner);
        Info << " goes from " << min(T).value() << " to " << max(T).value()
             << endl;
        T.write();

        Info << "Calculating density, rho" << flush;
        volScalarField rho("rho", p/(R*T));
        Info << " goes from " << min(rho).value() << " to " << max(rho).value()
             << endl;
        rho.write();
        
        Info << "Calculating Brunt Vaiasalla frequency" << flush;
        volScalarField BruntFreq
        (
            IOobject("BruntFreq", runTime.timeName(), mesh),
            Foam::sqrt(-(g & fvc::grad(theta))/theta)
        );
        Info << " goes from " << min(BruntFreq).value() << " to "
             << max(BruntFreq).value() << endl;
        BruntFreq.write();
    }

    Info << endl;

    return(0);
}


// ************************************************************************* //
