/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    setJetVorticity

Description
    Set the initial vorticity for a jet defined by a quintic polynomial

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("dictionary name (in system)");
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"
#   include "createTime.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word dictName = args.args()[1].c_str();
    Info<< "Reading initial conditions from" << dictName << endl;

    IOdictionary initDict
    (
        IOobject
        (
            dictName, mesh.time().system(), mesh, IOobject::MUST_READ
        )
    );

    // Maximum jet velocity
    const dimensionedScalar u0(initDict.lookup("Umax"));
    // y location of the centre of the jet
    const dimensionedScalar yc(initDict.lookup("jetCentre"));
    // jet half width
    const dimensionedScalar w(initDict.lookup("jetHalfWidth"));
    
    const scalar yb = (yc - w).value();
    const scalar yt = (yc + w).value();

    volScalarField y("y", mesh.C().component(vector::Y));
    volScalarField yy("yy", (y - yc)/w);

    // Create initial vorticity
    volScalarField vorticity
    (
        IOobject("vorticity", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    vorticity = 6*yy*u0/w*(pow(yy,4) - 2*sqr(yy) + 1);
    
    forAll(vorticity, i)
    {
        if (y[i] < yb || y[i] > yt)
        {
            vorticity[i] = 0;
        }
    }
    vorticity.correctBoundaryConditions();
    
    vorticity.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
