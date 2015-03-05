/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    setHorizontalVelocityField

Description
    Sets the horizontal velocity field for the 
    Schar et al Mon. Wea. Rev., 130(10):2459-2480, 2002
    horizontal advection over orography test case.
    The flux is set by integrating the stream function around each face
    to ensure that the flux is divergence free

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

class VelocityProfile
{
    public:
    VelocityProfile(const IOdictionary& dict)
        :
            u0(readScalar(dict.lookup("maxVelocity"))),
            z1(readScalar(dict.lookup("zeroVelocityHeight"))),
            z2(readScalar(dict.lookup("maxVelocityHeight")))
    {}

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyTo(GeometricField<Type, PatchField, GeoMesh>& field)
    {
        applyToInternalField(field);
        applyToBoundary("inlet", field);
        applyToBoundary("outlet", field);
        applyToBoundary("top", field);
        applyToBoundary("ground", field);
    }

    const scalar u0;
    const scalar z1;
    const scalar z2;
    
    private:
    void applyToInternalField(surfaceVectorField& field)
    {
        forAll(field, faceI)
        {
            const point& Cf = field.mesh().Cf()[faceI];
            field[faceI] = velocityAt(Cf.z());
        }
    }

    void applyToInternalField(volVectorField& field)
    {
        forAll(field, cellI)
        {
            const point& C = field.mesh().C()[cellI];
            field[cellI] = velocityAt(C.z());
        }
    }

    template<template<class> class PatchField, class GeoMesh>
    void applyToBoundary
    (
        const word name, GeometricField<vector, PatchField, GeoMesh>& field
    )
    {
        const label boundaryI = field.mesh().boundaryMesh().findPatchID(name);
        forAll(field.boundaryField()[boundaryI], cellI)
        {
            const point& face = field.mesh().Cf().boundaryField()[boundaryI][cellI];
            field.boundaryField()[boundaryI][cellI] = velocityAt(face.z());
        }
    }
    
    //- The flux field from the streamfunction by application of Stokes' theorem
    //- for the internal flux field
    void applyToInternalField(surfaceScalarField& phi)
    {
        const fvMesh& mesh = phi.mesh();
        forAll(phi, faceI)
        {
            // Circulate around the vertices of the face and sum to contribution
            // to the flux
            const face& f = mesh.faces()[faceI];
            point p0 = mesh.points()[f.last()];
            point p1 = mesh.points()[f.first()];
            point pmid = 0.5*(p0 + p1);
            phi[faceI] = streamFunctionAt(pmid.z())*(p1.y() - p0.y());
            for(label ip = 1; ip < f.size(); ip++)
            {
                p0 = p1;
                p1 = mesh.points()[f[ip]];
                point pmid = 0.5*(p0 + p1);
                phi[faceI] += streamFunctionAt(pmid.z())*(p1.y() - p0.y());
            }
        }
    }
    
    //- The flux field from the streamfunction by application of Stokes' theorem
    //- for the boundary
    void applyToBoundary(const word name, surfaceScalarField& phi)
    {
        const fvMesh& mesh = phi.mesh();
        const label patchI = mesh.boundaryMesh().findPatchID(name);
        scalarField& bf = phi.boundaryField()[patchI];
        forAll(bf, faceI)
        {
            const face& f = mesh.boundaryMesh()[patchI][faceI];
            point p0 = mesh.points()[f.last()];
            point p1 = mesh.points()[f.first()];
            point pmid = 0.5*(p0 + p1);
            bf[faceI] = streamFunctionAt(pmid.z())*(p1.y() - p0.y());
            for(label ip = 1; ip < f.size(); ip++)
            {
                p0 = p1;
                p1 = mesh.points()[f[ip]];
                point pmid = 0.5*(p0 + p1);
                bf[faceI] += streamFunctionAt(pmid.z())*(p1.y() - p0.y());
            }
        }
    }

    vector velocityAt(const scalar z) const
    {
        if (z > z1 && z < z2)
        {
            return vector(u0*sqr(Foam::sin(M_PI/2*(z-z1)/(z2-z1))), 0, 0);
        }
        else if (z >= z2)
        {
            return vector(u0, 0, 0);
        }
        else
        {
            return vector(0, 0, 0);
        }
    }
    
    scalar streamFunctionAt(const scalar z) const
    {
        if (z <= z1) return 0;
        else if (z <= z2)
        {
            return -0.5*u0*(z - z1 - (z2-z1)/M_PI*Foam::sin(M_PI*(z-z1)/(z2-z1)));
        }
        else return -0.5*u0*(2*z - z2 - z1);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

    Info << "Reading velocityFieldDict" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    VelocityProfile velocityProfile(initDict);
    Info << "Creating velocity field Uf" << endl;
    velocityProfile.applyTo(Uf);
    Uf.write();

    Info << "Creating flux field, phi" << endl;
    velocityProfile.applyTo(phi);
    phi.write();
    
    Info << "Calculating the divergence field to check that it is zero" << endl;
    volScalarField divu("divu", fvc::div(phi));
    divu.write();
    
    Info << "Correcting Uf based on the flux" << endl;
    Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    Uf.write();

    Info << "Creating velocity field U" << endl;
    velocityProfile.applyTo(U);
    U.write();
}
