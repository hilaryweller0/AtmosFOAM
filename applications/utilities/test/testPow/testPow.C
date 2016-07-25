/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    testPow

Description
    Test pow(volScalarField, volScalarField) and associated dimensions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Calculate a^b where a has dimensions mass and b = 2.0

    volScalarField a
    (
        IOobject("a", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("a", dimensionSet(1,0,0,0,0), scalar(1))
    );
    volScalarField b
    (
        IOobject("b", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("b", dimless, scalar(2))
    );
    volScalarField powab("powab", pow(a,b));
    volScalarField aSquared("aSquared", pow(a,2));
    
    Info << "a has dimensions " << a.dimensions() << nl
         << "b has dimensions " << b.dimensions() << nl
         << "a^b has dimensions " << powab.dimensions() << nl
         << "whereas a^2 has dimensions " << aSquared.dimensions() << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
