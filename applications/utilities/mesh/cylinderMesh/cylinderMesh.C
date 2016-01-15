/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    cylinderMesh

Description
    Takes a flat 2d mesh in the x-z plane and transforms it into a cylinder. 
    You then need to run stitchMesh to joint the boundary faces together. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The width of the mesh in the x direction (which will become the diameter
    // of the cylinder
    const boundBox bounds(mesh.points());
    const scalar width = bounds.max()[0] - bounds.min()[0];
    
    // the depth of the mesh (the y bounds) in order to set the radius of the
    // inside and outside points of the cylinder
    const scalar depth = bounds.max()[1] - bounds.min()[1];
    const scalar rInside  = 0.5*width - 0.5*depth;
    const scalar rOutside = 0.5*width + 0.5*depth;

    // Create a new set of points that will wrap the mesh around a cylinder
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );
    
    // Wrap all points around the cylinder
    const scalar small = SMALL*bounds.span()[1];
    const scalar meshymin = bounds.min()[1];
    const scalar meshymax = bounds.max()[1];
    forAll(newPoints, ip)
    {
        point& p = newPoints[ip];
        // Check to see if the point is on the front or back patch
        
        if (mag(p.y() - meshymax) < small)
        {
            p = point
            (
                rInside*Foam::cos(2*M_PI*p.x()/width),
                rInside*Foam::sin(2*M_PI*p.x()/width),
                M_PI*((p.z() - bounds.min()[2])/bounds.span()[2]-0.5)
            );
        }
        else if (mag(p.y() - meshymin) < small)
        {
            p = point
            (
                rOutside*Foam::cos(2*M_PI*p.x()/width),
                rOutside*Foam::sin(2*M_PI*p.x()/width),
                M_PI*((p.z() - bounds.min()[2])/bounds.span()[2]-0.5)
            );
        }
        else
        {
            FatalErrorIn("cylinderMesh")
                << " can only transform a mesh which is 2d in the xz plane.\n"
                << "This mesh goes on the y direction from "
                << bounds.max()[0] << " to " << bounds.max()[1]
                << "\nbut point " << ip << " is at location " << p
                << exit(FatalError);
        }
    }
    newPoints.write();

    Info<< "End\n" << endl;
}
