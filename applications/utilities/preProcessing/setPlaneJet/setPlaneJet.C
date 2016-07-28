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
    setPlaneJet

Description
    Set the initial U and h for a shallow water eqn jet on a beta plane

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

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
            "initialConditions",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Maximum jet velocity
    const dimensionedScalar u0(initDict.lookup("u0"));
    // y location of the centre of the jet
    const dimensionedScalar yc(initDict.lookup("jetCentre"));
    // jet half width
    const dimensionedScalar w(initDict.lookup("jetHalfWidth"));
    // Shallow water height at the equator
    const dimensionedScalar h0(initDict.lookup("h0"));
    
    // Height, location and widths of a bump
    const dimensionedScalar bumpHeight
    (
        initDict.lookupOrDefault<dimensionedScalar>
        (
            "bumpHeight", dimensionedScalar("bumpHeight", dimLength, scalar(0))
        )
    );
    const dimensionedVector bumpCentre
    (
        initDict.lookupOrDefault<dimensionedVector>
        (
            "bumpCentre", dimensionedVector("bumpCentre", dimLength, point::zero)
        )
    );
    const dimensionedVector bumpWidth
    (
        initDict.lookupOrDefault<dimensionedVector>
        (
            "bumpWidth", dimensionedVector("bumpWidth", dimLength, point::zero)
        )
    );
    
    const scalar yb = (yc - w).value();
    const scalar yt = (yc + w).value();

    Info<< "Reading environmental properties\n" << endl;

    IOdictionary envDict
    (
        IOobject
        (
            "environmentalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // gravity magnitude
    const dimensionedScalar g(envDict.lookup("g"));
    // beta
    const dimensionedScalar beta(envDict.lookup("beta"));

    // Create initial height and velocity field
    volScalarField h
    (
        IOobject("h", runTime.timeName(), mesh),
        mesh,
        h0,
        "zeroGradient"
    );
    
    volScalarField y("y", mesh.C().component(vector::Y));
    volScalarField yy("yy", (y - yc)/w);
    volScalarField u("u", u0*(1 - 3*sqr(yy) + 3*pow(yy,4) - pow(yy,6)));
    u = max(u, u0*0);

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        vector(1.,0.,0.)*u,
        "slip"
    );
        
    surfaceScalarField yf("yf", mesh.Cf().component(vector::Y));
    surfaceScalarField yyf("yyf", (yf - yc)/w);
    surfaceScalarField uf("uf", u0*(1 - 3*sqr(yyf) + 3*pow(yyf,4) - pow(yyf,6)));
    uf = max(uf, u0*0);

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        vector(1.,0.,0.)*uf
    );
    
    dimensionedScalar hmin = h0 - 32/35.*u0*beta/g*w*yc;
    
    forAll(h, i)
    {
        if (y[i] > yb && y[i] < yt)
        {
            dimensionedScalar hh = h0 - u0*beta/g*w*
            (
                -0.125*w + 16/35.*yc
              + yc*yy[i]*
                (1 - sqr(yy[i]) + 0.6*pow(yy[i],4) - 1/7.*pow(yy[i],6))
              + w*sqr(yy[i])*
                (0.5 - .75*sqr(yy[i]) + 0.5*pow(yy[i],4) - 0.125*pow(yy[i],6))
            );
        
            h[i] = hh.value();
        }
        else if(y[i] >= yt)
        {
            h[i] = hmin.value();
        }
    }
    h.correctBoundaryConditions();
    
    if (bumpHeight.value() > SMALL)
    {
        volScalarField hun("hunperturbed", h);
        hun.write();
        
        volScalarField x("x", mesh.C().component(vector::X));
        
        forAll(h, i)
        {
            if (y[i] > yb && y[i] < yt)
            {
                h[i] += bumpHeight.value()
                       *sqr(Foam::cos(0.5*pi*yy[i]))
                       *Foam::exp(-sqr((x[i]-bumpCentre.value()[0])/bumpWidth.value()[0]))
                       *Foam::exp(-sqr((y[i]-bumpCentre.value()[1])/bumpWidth.value()[1]));
            }
        }
    }
    
    h.write();
    U.write();
    Uf.write();
    volVectorField hU
    (
        IOobject("hU", runTime.timeName(), mesh),
        h*U,
        "slip"
    );
    hU.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
