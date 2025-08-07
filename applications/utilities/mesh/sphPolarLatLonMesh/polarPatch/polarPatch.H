// The FOAM Project // File: polarPatch.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   polarPatch
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    polarPatch

DESCRIPTION
    Defines a surface mesh with regular latitude-longitude vertices over a 
    sphere with un-refinements in the longitude direction so as to maintain a
    more constant spacing of meridians
    maxDxLonRatio = 0 means regular grid
    maxDxLonRatio = 2 means the number longitudes halves when the ratio 
                      of longitude distance to equatorial longitude
                      distance reaches a factor of 2

C++
    listOfSourceFiles
        polarPatch.C
    endListOfSourceFiles

AUTHOR
    Hilary Spencer.
*/
// ------------------------------------------------------------------------- //

#ifndef polarPatch_H
#define polarPatch_H

#include "PrimitivePatch.H"
#include "faceList.H"
#include "pointField.H"
#include "polarPatchData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polarPatch Declaration
\*---------------------------------------------------------------------------*/

class polarPatch
:
    public polarPatchData,
    public PrimitivePatch<faceList, pointField>
{
    // Private data

        const scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        polarPatch(const polarPatch&);

        //- Disallow default bitwise assignment
        void operator=(const polarPatch&);


    // Helper functions

        //- return the list of faces
        faceList calcPolarFaces(const polarPatchData& grid) const;

        //-  return the vertex locations
        pointField calcVerts(const polarPatchData& grid, const scalar r) const;

        //- correct the direction of the faces
        void correctDirection();


public:

    // Constructors

        //- Construct  for a regular latlon grid and from radius, maxDxLonRatio
        //  nlon being the maximum number of longitudes for any latitude.
        //  Includes a single polyhedral cell at each pole
        //  maxDxLonRatio = 0 means regular grid
        //  maxDxLonRatio = 2 means the number longitudes halves when the ratio 
        //                    of longitude distance to equatorial longitude
        //                    distance reaches a factor of 2
        polarPatch
        (
            const label nlat, 
            const label nlon,
            const scalar radius,
            const scalarList& urefLats,
            const scalar maxDxLonRatio,
            const vector& meshRotation,
            const bool polarCell
        );
        
        //- Contruct from dictionary
        polarPatch
        (
            const IOdictionary& earthProperties
        );

    // Member functions
    
        //- Calculate the Voronoi Points for each face
        vectorField voronoiPoints() const;

        //- Calculate the edge Points for each face given the Voronoi Points
        vectorField edgePoints(const vectorField& vPts) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
