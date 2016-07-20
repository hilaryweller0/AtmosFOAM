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
#include "moistThermo.H"
#include "OFstream.H"

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
    #include "readThermoProperties.H"
    #include "readThermoPropertiesMoist.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        //#   include "createMesh.H"
        mesh.readUpdate();

        Info << "Reading Exner function" << endl;
        volScalarField Exner
        (
            IOobject("Exner", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        Info << "Reading thetaRho" << endl;
        volScalarField thetaRho
        (
            IOobject("thetaRho", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        Info << "Reading qv" << endl;
        volScalarField qv
        (
            IOobject("qv", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        Info << "Reading ql" << endl;
        volScalarField ql
        (
            IOobject("ql", runTime.timeName(),mesh,IOobject::MUST_READ),
            mesh
        );

        // Related mixing ratios
        volScalarField qd("qd", 1 - qv - ql);
        volScalarField rv("rv", qv/qd);
        volScalarField rl("rl", ql/qd);
        rv.write();
        rl.write();

        // Thermodynamics variable kappa as a function of mixing ratios
        volScalarField cpml("cpml", Cp + Cpv*rv + rl*Cpl);
        volScalarField Rm("Rm", R + Rv*rv);
        volScalarField kappam("kappam", Rm/cpml);

        Info<< "Creating field rho from eqn of state\n" << endl;
        volScalarField rho
        (
            IOobject("rho", runTime.timeName(), mesh),
            moistGasExnerRho(Exner, thetaRho, kappam, R, pRef)
        );
        rho.write();

        Info << "Calculating pressure, p" << flush;
        volScalarField p("p", pFromExner(Exner, kappam, pRef));
        Info << " goes from " << min(p).value() << " to " << max(p).value()
             << endl;
        p.write();

        Info << "Calculating temperature, T" << flush;
        volScalarField T("T", TFromThetaRho(thetaRho, Exner, rv, rl, epsilon));
        Info << " goes from " << min(T).value() << " to " << max(T).value()
             << endl;
        T.write();

        Info << "Calculating thetae" << flush;
        volScalarField Lv("Lv", latentHeat(T, Lv0, T0, Cpl, Cpv));
        volScalarField thetae
        (
            IOobject
            (
                "thetae", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::AUTO_WRITE
            ),
            thetaeFromPrimitive(T, p, Lv, rv, rl, pRef, epsilon, Cp, Cpl, R)
        );
        Info << " goes from " << min(thetae).value() << " to "
             << max(thetae).value() << endl;
        thetae.write();
    }

    Info << endl;

    return(0);
}


// ************************************************************************* //
