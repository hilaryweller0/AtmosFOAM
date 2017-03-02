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
    testFVFCFstencil

Description
    Finds faces that share a vertex with another face

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace fv;

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        const face& f = mesh.faces()[faceI];
        labelList connectedFaceList;

        forAll(f, pointForFaceI)
        {
            const label pointI = f[pointForFaceI];
            const labelList& connectedFaces = mesh.pointFaces()[pointI];

            forAll(connectedFaces, i)
            {
                connectedFaceList.append(connectedFaces[i]);
            }
        }
        HashSet<label, Hash<label>> connectedFaceSet(connectedFaceList);

        Info << "faces connected to faceI " << faceI << " are " << connectedFaceSet << endl;
    }
}
