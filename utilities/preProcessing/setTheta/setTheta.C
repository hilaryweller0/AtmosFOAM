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
    setTheta

Description
    Set theta based on an array of Brunt Vailsalla frequencies for different
    layers

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    // The Brunt Vaiasalla freqencies for the different layers
    const scalarList Brunt(envProperties.lookup("BruntVaisallaFreq"));
    // The extents of the different layers (size one greater)
    const scalarList zN(envProperties.lookup("zN"));
    if (Brunt.size()+1 != zN.size())
    {
        FatalErrorIn("setTheta")
            << " size of BruntVaisallaFreq in environmentalProperties should be"
            << " one smaller than the size of zN"
            << exit(FatalError);
    }
        
    Info<< "Reading theta_init\n" << endl;
    volScalarField theta_init
    (
        IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info<< "Creating theta\n" << endl;
    volScalarField theta
    (
        IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
        theta_init
    );
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Loop over the different layers and set theta
    forAll(theta, cellI)
    {
        const scalar z = mesh.C()[cellI].z();
        bool found = false;
        scalar theta0 = T0.value();
        for(label il = 0; il < Brunt.size(); il++)
        {
            if (z >= zN[il] && z < zN[il+1])
            {
                theta[cellI] = theta0*Foam::exp
                (
                    sqr(Brunt[il])/mag(g.value())*(z-zN[il])
                );
                found = true;
            }
            theta0 *= Foam::exp
                (sqr(Brunt[il])/mag(g.value())*(zN[il+1]-zN[il]));
        }
        if (!found)
        {
            FatalErrorIn("setTheta") << " cell " << cellI
                << " with height " << z << " not found in levels "
                << zN << exit(FatalError);
        }
    }
    forAll(theta.boundaryField(), patchI)
    {
        fvPatchField<scalar>& thetap = theta.boundaryField()[patchI];
        forAll(thetap, facei)
        {
            const scalar z = mesh.C().boundaryField()[patchI][facei].z();
            bool found = false;
            scalar theta0 = T0.value();
            for(label il = 0; il < Brunt.size(); il++)
            {
                if (z >= zN[il] && z < zN[il+1])
                {
                    thetap[facei] = theta0*Foam::exp
                    (
                        sqr(Brunt[il])/mag(g.value())*(z-zN[il])
                    );
                    found = true;
                }
                theta0 *= Foam::exp
                    (sqr(Brunt[il])/mag(g.value())*(zN[il+1]-zN[il]));
            }
            if (!found)
            {
                FatalErrorIn("setTheta") << " boundary face " << facei
                    << " with height " << z << " not found in levels "
                    << zN << exit(FatalError);
            }                
        }
        
    }

    theta.write();

    surfaceScalarField thetaf("thetaf", fvc::interpolate(theta));
    thetaf.write();

    volScalarField BruntFreq
    (
        IOobject("BruntFreq", runTime.timeName(), mesh),
        Foam::sqrt(-(g & fvc::grad(thetaf))/theta)
    );
    BruntFreq.write();

    surfaceScalarField BruntFreq2f
    (
        IOobject("BruntFreq2f", runTime.timeName(), mesh),
        -(g & mesh.delta())/mag(mesh.delta())*fvc::snGrad(theta)/thetaf
    );
    BruntFreq2f.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

