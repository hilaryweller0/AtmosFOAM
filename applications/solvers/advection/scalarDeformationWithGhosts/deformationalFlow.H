/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::deformationalFlow

Description
    Reads in the parameters that define the deformational flow and updates
    the flux and velocity field in member functions

\*---------------------------------------------------------------------------*/

#ifndef deformationalFlow_H
#define deformationalFlow_H

#include "wallFvPatch.H"
#include "coupledFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class deformationalFlow Declaration
\*---------------------------------------------------------------------------*/

class deformationalFlow
{
    // Private data

        //- Constants in the specification of the velocity field
        const scalar endTime_;
        const scalar psiScale_;
        const vector L_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        deformationalFlow(const deformationalFlow&);

        //- Disallow default bitwise assignment
        void operator=(const deformationalFlow&);


public:

    // Constructors

        //- Construct from dictionary
        deformationalFlow(IOdictionary initDict)
        :
            endTime_(dimensionedScalar(initDict.lookup("endTime")).value()),
            psiScale_(readScalar(initDict.lookup("psiScale"))),
            L_(dimensionedVector(initDict.lookup("L")).value())
        {}

    //- Destructor
    ~deformationalFlow()
    {}


    // Member Functions

    //- Update the flux based on the current time
    void update(surfaceScalarField& phi)
    {
        // Get the current time and the mesh
        const fvMesh& mesh = phi.mesh();
        const scalar t = mesh.time().value();

        // Set internal values of phi
        phi.ref() = 0;
        forAll(phi, faceI)
        {
            const face& f = mesh.faces()[faceI];
            point p0 = mesh.points()[f.last()];
            point p1 = p0;
            
            // Loop around all edges of the face to set the flux from the 
            // streamfunction using Stokes circulation theorem
            forAll(f, ip)
            {
                p0 = p1;
                p1 = mesh.points()[f[ip]];
                point pmid = 0.5*(p0 + p1);
                scalar z = pmid.z();
                scalar x = pmid.x();
                vector streamFunction = vector(0,-1,0)*
                (
                    psiScale_*sqr(0.5*L_.x()/M_PI)/endTime_*sqr
                    (
                        sin(2*M_PI*(x/L_.x() - t/endTime_))
                       *cos(M_PI*z/L_.z())
                    )
                    *cos(M_PI*t/endTime_)
                  - L_.x()*z/endTime_
                );
                phi[faceI] += streamFunction & (p0 - p1);
            }
        }
        
        // Set boundary values of phi (if needed)
        forAll(phi.boundaryField(), patchI)
        {
            if
            (
                phi.boundaryField()[patchI].type() == "calculated"
             || phi.boundaryField()[patchI].type() == "cyclic"
            )
            {
                forAll(phi.boundaryField()[patchI], fi)
                {
                    // Initialise the flux to zero before adding contributions
                    phi.boundaryFieldRef()[patchI][fi] = 0;
                
                    // Loop around edges of the face and set the flux
                    const face& f = mesh.boundaryMesh()[patchI][fi];
                    point p0 = mesh.points()[f.last()];
                    point p1 = p0;
                    forAll(f, ip)
                    {
                        p0 = p1;
                        p1 = mesh.points()[f[ip]];
                        point pmid = 0.5*(p0 + p1);
                        scalar z = pmid.z();
                        scalar x = pmid.x();
                        vector streamFunction = vector(0,-1,0)*
                        (
                            psiScale_*sqr(0.5*L_.x()/M_PI)/endTime_*sqr
                            (
                                sin(2*M_PI*(x/L_.x() - t/endTime_))
                               *cos(M_PI*z/L_.z())
                            )
                            *cos(M_PI*t/endTime_)
                          - L_.x()*z/endTime_
                        );
                        phi.boundaryFieldRef()[patchI][fi] 
                            += streamFunction & (p0 - p1);
                    }
                }
            }
        }
    }

    //- Update the flux and velocity based on the current time
    void update
    (
        surfaceScalarField& phi,
        volVectorField& U,
        surfaceVectorField& Uf
    )
    {
        const fvMesh& mesh = phi.mesh();
        update(phi);
        U = fvc::reconstruct(phi);
        Uf = linearInterpolate(U);
        Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
