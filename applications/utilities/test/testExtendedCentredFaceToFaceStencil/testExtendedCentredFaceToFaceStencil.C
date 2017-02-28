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
    testExtendedCentredFaceToFaceStencil

Description
    Applies the stencil to a mesh to check what the stencil looks like.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "extendedCentredFaceToFaceStencil.H"
#include "CFCFaceToFaceStencil.H"

using namespace fv;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    surfaceScalarField Tf
    (
        IOobject
        (
            "Tf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    CFCFaceToFaceStencil f2fStencil(mesh);
    extendedCentredFaceToFaceStencil stencil(f2fStencil);

    Info << "stencil faceI\n" << stencil.elements() << endl;

    List<List<point> > stencilPoints(mesh.nFaces());

    stencil.collectData
    (
        mesh.Cf(),
        stencilPoints
    );

    Info << "stencil points\n" << stencilPoints << endl;

    List<List<vector > > coeffs(mesh.nFaces());

    forAll(stencilPoints, stencilForFaceI)
    {
        forAll(stencilPoints[stencilForFaceI], faceIInStencil)
        {
            coeffs[stencilForFaceI].append(vector(1, 0, 0));
        }
    }

    tmp<surfaceVectorField> tsvf = stencil.weightedSum(Tf, coeffs);
    Info << "grad(Tf)" << tsvf() << endl;

    // for each stencil, we want to calculate the gradient which is a weighted sum
    // can we use stencil.weightedSum()?  yes! I think so.
}


// ************************************************************************* //
