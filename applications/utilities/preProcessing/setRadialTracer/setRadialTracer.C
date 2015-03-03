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
    setRadialTracer

Description
    Set an elliptical tracer, T, on a 2D x-z mesh.  Tracer magnitude is at its
    maximum at the centre and zero at its perimeter.
    
    Expects the file $CASE/constant/radialTracerDict to contain:
    - centre_x, centre_z          The (x,z) centre of the tracer
    - halfWidth_x, halfWidth_z    The tracer half widths in the x and z coords
    - max_magnitude               The maximum tracer magnitude at its centre

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    Foam::argList::addOption("x0", "int", "specify horizontal placement of tracer, overrides x0 specified in radialTracerDict");
    Foam::argList::addOption("fieldName", "filename", "specify the name of the tracer field name (default 'T')");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info << "Reading radialTracerDict" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "radialTracerDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const scalar rho0(readScalar(dict.lookup("max_magnitude")));
    scalar x0(readScalar(dict.lookup("centre_x")));
    const scalar z0(readScalar(dict.lookup("centre_z")));
    const scalar Ax(readScalar(dict.lookup("halfWidth_x")));
    const scalar Az(readScalar(dict.lookup("halfWidth_z")));

    string fieldName = "T";
    if (args.options().found("fieldName")) {
        fieldName = args.options()["fieldName"];
    }

    if (args.options().found("x0")) {
        x0 = readScalar(IStringStream(args.options()["x0"])());
    }

    Info << "Creating tracer field " << fieldName << endl;

    volScalarField T
    (
        IOobject(fieldName, runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(fieldName, dimless, scalar(0)),
        "zeroGradient"
    );

    forAll(T, cellI) {
        const point& c = mesh.C()[cellI];
        
        double r = Foam::sqrt(pow((c.x()-x0)/Ax,2)+pow((c.z()-z0)/Az,2));
        
        if (r<=1) {
            T[cellI] = rho0*pow(Foam::cos(M_PI*r/2),2);
        }       
    };
    
    T.write();
    
    Info<< "End\n" << endl;
    
    return(0);
}
