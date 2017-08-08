// The FOAM Project // File: sphPolarMesh.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   sphPolarMesh
   \\  /           | Family: util
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR
    Henry G. Weller.

-------------------------------------------------------------------------------
*/

#include "Time.H"

#include "argList.H"
#include "polarPatch.H"
#include "dimensionedTypes.H"
#include "mathematicalConstants.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
#include "OFstream.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "readEarthProperties.H"

    polarPatch ppatch(earthProperties);

    extrudedMesh eMesh
    (
        IOobject
        (
            "",
            earthProperties.time().constant(),
            earthProperties.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
//        runTime,
        ppatch,
        extrudeModels::linearRadial(earthProperties)
    );

    if (!eMesh.write())
    {
        FatalErrorIn(args.executable()) << "Failed writing mesh"
            << exit(FatalError);
    }

    //eMesh.checkMesh();
    
//    vectorField VoronoiPoints(ppatch.voronoiPoints());
//    OFstream osv("VoronoiPoints");
//    osv << VoronoiPoints;
//    
//    vectorField edgePoints(ppatch.edgePoints(VoronoiPoints));
//    OFstream ose("edgePoints");
//    ose << edgePoints;
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
