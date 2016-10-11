/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    writeSphericalMeshGeometry

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Total mesh volume = " << sum(mesh.V()) << endl;

    const fileName outFile = args.rootPath() / args.caseName() / runTime.constant();

    const fileName pointFile = outFile / "point.diagnostics";
    OFstream pos(pointFile);
    pos << "# pointI x y z mag(position)" << endl;
    forAll(mesh.points(), pointI)
    {
        pos << pointI << " " << mesh.points()[pointI][0] << " "
            << mesh.points()[pointI][1] << " "
            << mesh.points()[pointI][2] << " "
            << mag(mesh.points()[pointI]) << endl;
    }

    const fileName cellFile = outFile / "cell.diagnostics";
    OFstream cos(cellFile);
    cos << "# cellI x y z mag(C) volume" << endl;
    for(label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        cos << cellI << " " << mesh.C()[cellI][0] << " " << mesh.C()[cellI][1] << " "
            << mesh.C()[cellI][2] << " " << mag(mesh.C()[cellI]) << " "
            << mesh.V()[cellI] << endl;
    }

    const fileName faceFile = outFile / "face.diagnostics";
    OFstream fos(faceFile);
    fos << "# faceI x y z mag(Cf) area radialComponent kComponent" << endl;
    
    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        scalar radialComponent = (mesh.Cf()[faceI]/mag(mesh.Cf()[faceI])) & (mesh.Sf()[faceI]/mesh.magSf()[faceI]);
        scalar kComponent = mesh.Sf()[faceI] & vector(0,0,1);
        
        fos << faceI << " " << mesh.Cf()[faceI][0] << " " 
           << mesh.Cf()[faceI][1] << " "
           << mesh.Cf()[faceI][2] << " " << mag(mesh.Cf()[faceI]) << " "
           << mesh.magSf()[faceI] << " " << radialComponent << " "
           << kComponent << endl;
    }

    return EXIT_SUCCESS;
}
