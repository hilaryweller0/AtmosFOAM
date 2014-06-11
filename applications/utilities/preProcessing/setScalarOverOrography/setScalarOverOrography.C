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
    setScalarOverOrography

Description
    Set the velocity, U, and the initial scalar, T for scalar transport over
    orography

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addOption("x0", "int", "specify horizontal placement of tracer, overrides x0 specified in initialConditions dictionary");
    Foam::argList::addBoolOption("withoutWindField", "omit wind field U from output");
    Foam::argList::addOption("tracerFieldFileName", "filename", "specify the name of the tracer field file name (default 'T')");
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "Reading initial conditions" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "initialConditions",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Maximum wind speed
    const scalar u0(readScalar(initDict.lookup("u0")));
    // Height at which the winds attains its maximum value
    const scalar z2(readScalar(initDict.lookup("z2")));
    // Height below which the winds are zero
    const scalar z1(readScalar(initDict.lookup("z1")));
    // Initial maximum tracer value
    const scalar rho0(readScalar(initDict.lookup("rho0")));
    // Initial tracer position
    scalar x0(readScalar(initDict.lookup("x0")));
    const scalar z0(readScalar(initDict.lookup("z0")));
    // Half widths
    const scalar Ax(readScalar(initDict.lookup("Ax")));
    const scalar Az(readScalar(initDict.lookup("Az")));
    string tracerFieldFileName = "T";
    if (args.options().found("tracerFieldFileName")) {
        tracerFieldFileName = args.options()["tracerFieldFileName"];
    }

    if (args.options().found("x0")) {
        x0 = readScalar(IStringStream(args.options()["x0"])());
    }
    
    Info << "Creating initial tracer field " << tracerFieldFileName << endl;
    volScalarField T
    (
        IOobject(tracerFieldFileName, runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(tracerFieldFileName, dimless, scalar(0)),
        "zeroGradient"
    );
    
    if (!args.options().found("withoutWindField")) {
	Info << "Creating initial wind field U" << endl;
	volVectorField U
	(
		IOobject("U", runTime.timeName(), mesh),
		mesh,
		dimensionedVector("U", dimVelocity, vector(0,0,0)),
		"zeroGradient"
	);

        forAll(T, cellI) {
            const point& c = mesh.C()[cellI]; // Gets the mesh values from mesh.C

	    if (c.z() > z1 && c.z() < z2) { // region of changing wind speed
		    U[cellI] = vector( pow((Foam::sin(M_PI/2*(c.z()-z1)/(z2-z1))),2), 0, 0 );
	    } else if (c.z() >= z2) { // region of constant max wind speed
		    U[cellI] = vector(u0, 0, 0);
	    }
	}

    	U.write();
    }
        
    // Loop through all T and set correct values depending on location
    forAll(T, cellI) {
        const point& c = mesh.C()[cellI]; // Gets the mesh values from mesh.C
        
        // Define r as used in the initial tracer field
        double r = Foam::sqrt(pow((c.x()-x0)/Ax,2)+pow((c.z()-z0)/Az,2));
        
        if (r<=1) {
            T[cellI] = rho0*pow(Foam::cos(M_PI*r/2),2);
        }       
    };
    
    T.write();
    
    Info<< "End\n" << endl;
    
    return(0);
}

// ************************************************************************* //
