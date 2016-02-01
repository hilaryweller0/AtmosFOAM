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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    testInterp

Description
    Tests the interpolation scheme on some polynomials. When testing an upwind
    interpolation scheme the flux should be defined as xf (or yf or zf). Eg:
    interpolationSchemes
    {
        default        cubicUpwindCPCFit xf 3;
    }

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace fv;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // A scale factor for the geometry
    dimensionedScalar a("length", dimLength, mesh.bounds().mag());

    // Non-dimensionalised coodinate directions
    volScalarField x
    (
        IOobject("x", runTime.timeName(), mesh),
        mesh.C().component(vector::X)/a
    );
    volScalarField y
    (
        IOobject("y", runTime.timeName(), mesh),
        mesh.C().component(vector::Y)/a
    );
    volScalarField z
    (
        IOobject("z", runTime.timeName(), mesh),
        mesh.C().component(vector::Z)/a
    );

    surfaceScalarField xf
    (
        IOobject("xf", runTime.timeName(), mesh),
        mesh.Cf().component(vector::X)/a
    );
    surfaceScalarField yf
    (
        IOobject("yf", runTime.timeName(), mesh),
        mesh.Cf().component(vector::Y)/a
    );
    surfaceScalarField zf
    (
        IOobject("zf", runTime.timeName(), mesh),
        mesh.Cf().component(vector::Z)/a
    );

    volScalarField triQuadratic = (x - y)*(x - z);
    triQuadratic.rename("triQuadratic");
    triQuadratic.write();

    volScalarField xSqr = sqr(x);
    xSqr.rename("xSqr");
    xSqr.write();

    volScalarField triCubic = x*triQuadratic;
    triCubic.rename("triCubic");
    triCubic.write();

    volScalarField xCube = xSqr*x;
    xCube.rename("xCube");
    xCube.write();

    surfaceScalarField triQuadf = fvc::interpolate(triQuadratic);
    triQuadf.rename("triQuadf");
    triQuadf.write();

    surfaceScalarField triCubf = fvc::interpolate(triCubic);
    triCubic.rename("triCubic");
    triCubic.write();

    surfaceScalarField xSqrf = fvc::interpolate(xSqr);
    xSqrf.rename("xSqrf");
    xSqrf.write();

    surfaceScalarField xCubef = fvc::interpolate(xCube);
    xCubef.rename("xCubef");
    xCubef.write();

    surfaceScalarField triQuadfa = (xf - yf)*(xf - zf);
    triQuadfa.rename("triQuadfa");
    triQuadfa.write();

    surfaceScalarField triCubfa = xf*triQuadfa;
    triCubfa.rename("triCubfa");
    triCubfa.write();

    surfaceScalarField xSqrfa = sqr(xf);
    xSqrfa.rename("xSqrfa");
    xSqrfa.write();

    surfaceScalarField xCubefa = xSqrfa*xf;;
    xCubefa.rename("xCubefa");
    xCubefa.write();

    surfaceScalarField triQuadError = triQuadf - triQuadfa;
    triQuadError.rename("triQuadError");
    triQuadError.write();

    surfaceScalarField triCubError = triCubf - triCubfa;
    triCubError.rename("triCubError");
    triCubError.write();

    surfaceScalarField xSqrError = xSqrf - xSqrfa;
    xSqrError.rename("xSqrError");
    xSqrError.write();

    surfaceScalarField xCubeError = xCubef - xCubefa;
    xCubeError.rename("xCubeError");
    xCubeError.write();

    scalarList magE(mag(triQuadError.internalField()));
    label maxIn(findMax(magE));
    Info << "maximum quadratic error " << triQuadError[maxIn]
        << " at face " << maxIn << " at location " << mesh.Cf()[maxIn] << endl;

    magE = scalarList(mag(triCubError.internalField()));
    maxIn = findMax(magE);
    Info << "maximum cubic error " << triCubError[maxIn]
        << " at face " << maxIn << " at location " << mesh.Cf()[maxIn] << endl;

    magE = scalarList(mag(xSqrError.internalField()));
    maxIn = findMax(magE);
    Info << "maximum x squared error " << xSqrError[maxIn]
        << " at face " << maxIn << " at location " << mesh.Cf()[maxIn] << endl;

    magE = scalarList(mag(xCubeError.internalField()));
    maxIn = findMax(magE);
    Info << "maximum x cubed error " << xCubeError[maxIn]
        << " at face " << maxIn << " at location " << mesh.Cf()[maxIn] << endl;
}


// ************************************************************************* //
