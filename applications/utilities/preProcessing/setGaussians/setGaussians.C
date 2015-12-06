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
    setGaussians.C

Description
    Initial a given volScalarField to be the sum of a set of Gaussians with 
    different maxima, radius and centres

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gaussian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading initial conditions\n" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "setGaussiansDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Read in background value
    dimensionedScalar value
    (
        initDict.lookupOrDefault<dimensionedScalar>
        (
            "backgroundValue",
            dimensionedScalar("T", dimless, scalar(0))
        )
    );

    // Read in list of Gaussians
    List<gaussian> gaussians(initDict.lookup("gaussians"));

    // Initialise tracer to the background value
    volScalarField T
    (
        IOobject
        (
            value.name(), runTime.timeName(), mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        value
    );
    
    // Add fields for all of the Gaussian distributions
    forAll(gaussians, ig)
    {
        T += gaussians[ig].field(mesh);
    }
    
    T.write();
    
    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
