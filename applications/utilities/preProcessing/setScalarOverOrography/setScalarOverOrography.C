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
    This utility is DEPRECATED.
    To set a horizontal velocity field, use setHorizontalVelocityField.
    To set a radial tracer field, use setRadialTracer.

    Set the velocity, U, and the initial scalar, T for scalar transport over
    orography

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

using namespace Foam::constant::mathematical;

scalar ScharCos(const scalar x, const scalar a) {
    return sqr(Foam::cos(M_PI*x/a));
}

scalar ScharCosSmooth(const scalar x, const scalar a, const scalar hm) {
    scalar h = 0;
    if (mag(x) < a)
    {
        h = hm*sqr(Foam::cos(0.5*M_PI*x/a));
    }
    return h;
}

Foam::scalar ScharExp(const scalar x, const scalar a, const scalar hm)
{
    return hm*Foam::exp(-sqr(x/a));
}

scalar height_schaerCos(scalar x, scalar a, scalar hm, scalar lambda) {
    return ScharCosSmooth(x, a, hm) * ScharCos(x, lambda);
}

scalar height_schaerExp(scalar x, scalar a, scalar hm, scalar lambda) {
    return ScharExp(x, a, hm) * ScharCos(x, lambda);
}

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
    Info << "This utility is DEPRECATED." << endl;
    Info << "To set a horizontal velocity field, use setHorizontalVelocityField." << endl;
    Info << "To set a radial tracer field, use setRadialTracer." << endl << endl;

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

    enum windFieldType{LINEAR, BTF};
    const word windFieldName(initDict.lookupOrDefault<word>("windFieldType", "LINEAR"));
    const windFieldType windField = windFieldName == "BTF" ? BTF : LINEAR;

    enum mountainType{SCHAER_COS, SCHAER_EXP};
    const word mountainTypeName(initDict.lookupOrDefault<word>("mountainType", "SCHAER_COS"));
    const mountainType mountain = mountainTypeName == "SCHAER_EXP" ? SCHAER_EXP : SCHAER_COS;
    
    Info << "Creating initial tracer field " << tracerFieldFileName << endl;
    volScalarField T
    (
        IOobject(tracerFieldFileName, runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(tracerFieldFileName, dimless, scalar(0)),
        "zeroGradient"
    );
    
    if (!args.options().found("withoutWindField")) {
        volVectorField U
        (
            IOobject("U", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("U", dimVelocity, vector(0,0,0)),
            "zeroGradient"
        );

        if (windField == LINEAR) {
            Info << "Creating initial linear wind field U" << endl;
            forAll(T, cellI) {
                const point& c = mesh.C()[cellI]; // Gets the mesh values from mesh.C

                if (c.z() > z1 && c.z() < z2) { // region of changing wind speed
                    U[cellI] = vector(u0*pow((Foam::sin(M_PI/2*(c.z()-z1)/(z2-z1))),2), 0, 0 );
                } else if (c.z() >= z2) { // region of constant max wind speed
                    U[cellI] = vector(u0, 0, 0);
                }
            }
        } else {
            const scalar zt(readScalar(initDict.lookup("zt")));
            const scalar a(readScalar(initDict.lookup("a")));
            const scalar hm(readScalar(initDict.lookup("hm")));
            const scalar lambda(initDict.lookupOrDefault<scalar>("lambda", scalar(0)));
            Info << "Creating initial BTF wind field U" << endl;

            if (mountain == SCHAER_COS) {
		    forAll(T, cellI) {
			const point& c = mesh.C()[cellI];

			scalar x = c.x();
			scalar z = c.z();
			scalar h = height_schaerCos(x, a, hm, lambda);
			scalar u = zt / (zt - h);

			scalar dhdx = - hm * pi * (1/(2*a)*pow(Foam::cos(pi * x/lambda), 2) * Foam::sin(pi*x/a) +
				pow(Foam::cos(pi * x / (2.0*a)), 2) * Foam::sin(2.0*pi*x/lambda)/lambda);

			if (x < -a || x > a) {
			    dhdx = 0.0;
			}
			scalar w = zt * dhdx * (zt - z) / pow(zt - h, 2);

			if (c.z() >= h) {
			    U[cellI] = vector(u*u0, 0, w*u0);
			}
		    }
            } else {
		    Info << zt << ' ' << a << ' ' << hm << ' ' << lambda << endl;
		    forAll(T, cellI) {
			const point& c = mesh.C()[cellI];

			scalar x = c.x();
			scalar z = c.z();
			scalar h = height_schaerExp(x, a, hm, lambda);
			scalar u = zt / (zt - h);

			scalar dhdx = - 2.0 * hm * Foam::exp(-pow(x/a, 2)) * Foam::cos(pi * x / lambda) * (
				pi * Foam::sin(pi*x/lambda)/lambda + x * Foam::cos(pi * x / lambda) / pow(a,2));

			scalar w = zt * dhdx * (zt - z) / pow(zt - h, 2);

			if (c.z() >= h) {
			    U[cellI] = vector(u*u0, 0, w*u0);
			}
		    }
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
