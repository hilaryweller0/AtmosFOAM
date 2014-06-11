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
    createSpongeLayer

Description
    Create a damping layer for verticaly pointing faces near the top of an 
    atmospheric model (optionaly also sponge at the inlet)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    Info << "\nReading environmentalProperties" << endl;

    IOdictionary envProperties
    (
        IOobject
        (
            "environmentalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    dimensionedVector g(envProperties.lookup("g"));
    dimensionedVector ghat = g/mag(g);

    Info << "Reading in sponge layer coefficients\n" << endl;
    const scalar zB(readScalar(envProperties.lookup("spongeBase")));
    const scalar zt(readScalar(envProperties.lookup("spongeTop")));
    const scalar muBar(readScalar(envProperties.lookup("spongeMean")));
    const scalar xSpongeCentre
    (
        envProperties.lookupOrDefault<scalar>("xSpongeCentre", scalar(0))
    );
    const scalar xSpongeEnd
    (
        envProperties.lookupOrDefault<scalar>("xSpongeEnd", scalar(0))
    );
    const scalar xSpongeLength = mag(xSpongeCentre - xSpongeEnd);
        
    Info<< "Creating muSponge\n" << endl;
    surfaceScalarField muSponge
    (
        IOobject("muSponge", runTime.constant(), mesh),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Loop over all faces and set muSponge
    forAll(muSponge, faceI)
    {
        // First check if face has a vertical normal
        if(mag(mesh.Sf()[faceI] ^ ghat.value()) < mesh.magSf()[faceI]*1e-6)
        {
            // height of face centre
            const scalar z = -(mesh.Cf()[faceI] & ghat.value());
            // x distance to x sponge centre
            const scalar xDist = mag(mesh.Cf()[faceI].x() - xSpongeCentre);
            
            // set the sponge value if the height is above sponge base
            if (z > zB)
            {
                muSponge[faceI] = muBar*sqr(Foam::sin(0.5*pi*(z-zB)/(zt-zB)));
            }
            else if (z > zt)
            {
                FatalErrorIn("createSpongeLayer") << "face " << faceI
                    << " has height " << z
                    << " but the sponge is defined to lie between " << zB
                    << " and " << zt << exit(FatalError);
            }
            
            // set the sponge value if x is between xMin and xMax
            if (xDist <= xSpongeLength)
            {
                muSponge[faceI] += muBar
                    *sqr(Foam::sin(0.5*pi*(xSpongeLength-xDist)/xSpongeLength));
            }
        }
    }   

    muSponge.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

