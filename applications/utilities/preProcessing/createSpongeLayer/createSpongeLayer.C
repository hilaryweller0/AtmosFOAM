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
    atmospheric model. muSponge is a surfaceScalarField.
    This is generalised to all boundaries - damping can be imposed for faces
    parallel to any boundary. 
    For co-located models, a diagonal volTensorField is created, muSpongeC.
    All boundaries are assumed to be in the x, y or z planes

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
    
    Info<< "Creating muSponge\n" << endl;
    surfaceScalarField muSponge
    (
        IOobject("muSponge", runTime.constant(), mesh),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    );
    volTensorField muSpongeC
    (
        IOobject("muSponge", runTime.constant(), mesh),
        mesh,
        dimensionedTensor("muSponge", dimless, tensor::zero)
    );

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

    const bool spongeOnCellCentres
    (
        readBool(envProperties.lookup("spongeOnCellCentres"))
    );

    const label nSponges(readLabel(envProperties.lookup("nSponges")));

    // Read in and set sponge coefficients for each sponge
    for(label iSponge = 0; iSponge < nSponges; iSponge++)
    {
        Info << "Reading in coefficients for sponge " << iSponge << endl;
        
        const vector spongeCentre
        (
            envProperties.lookup(word("spongeCentre"+std::to_string(iSponge)))
        );

        const vector spongeLength
        (
            envProperties.lookup(word("spongeLength"+std::to_string(iSponge)))
        );
        const vector spongeDir = spongeLength/(mag(spongeLength)+VSMALL);
        
        const scalar spongeMax(readScalar
        (
            envProperties.lookup(word("spongeMax"+std::to_string(iSponge)))
        ));
        

        // Add coefficients for this sponge for each sponge type

        // Loop over cells and set muSpongeC
        if (spongeOnCellCentres) forAll(muSpongeC, celli)
        {
            // Vector to sponge centre
            vector d = mesh.C()[celli] - spongeCentre;
            
            // Add sponge coeff if point close to spongeCentre in direction
            scalar dist = mag(d & spongeLength)/(spongeLength & spongeLength);
            if (dist < 1)
            {
                vector mu = spongeDir*spongeMax*sqr(Foam::cos(0.5*pi*dist));
                muSpongeC[celli] += tensor(mu.x(),0,0,0,mu.y(),0,0,0,mu.z());
            }
        }
        // Loop over all faces and set muSponge
        else forAll(muSponge, faceI)
        {
            // First check if face has a normal in direction spongeLength
            if(mag(mesh.Sf()[faceI] ^ spongeDir) < mesh.magSf()[faceI]*1e-6)
            {
                // Vector to sponge centre
                vector d = mesh.Cf()[faceI] - spongeCentre;
            
                // Add sponge coeff if point close to spongeCentre in direction
                scalar dist = mag(d & spongeLength)/(spongeLength & spongeLength);
                if (dist < 1)
                {
                    muSponge[faceI] += spongeMax*sqr(Foam::cos(0.5*pi*dist));
                }
            }
        }
    }

    if (spongeOnCellCentres) muSpongeC.write();
    else muSponge.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

