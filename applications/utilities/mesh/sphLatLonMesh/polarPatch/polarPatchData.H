// The FOAM Project // File: polarPatchData.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   polarPatchData
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    polarPatchData

DESCRIPTION
    Defines the regular latitude-longitude locations over a 
    sphere with un-refinements in the longitude direction so as to maintain a
    more constant spacing of meridians

C++
    listOfSourceFiles
        polarPatchData.C
    endListOfSourceFiles

AUTHOR
    Hilary Spencer.
*/
// ------------------------------------------------------------------------- //

#ifndef polarPatchData_H
#define polarPatchData_H

#include "scalarList.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

class IOdictionary;

/*---------------------------------------------------------------------------*\
                           Class polarPatchData Declaration
\*---------------------------------------------------------------------------*/

class polarPatchData
{
    // Private data

        const bool polarCell_;
        const scalar maxDxLonRatio_;
        scalarList lats_;
        const scalarList latCs_;
        const label jEq_;
        List<scalarList> lons_;
        const List<scalarList> lonCs_;
        const vector meshRotation_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        polarPatchData(const polarPatchData&);

        //- Disallow default bitwise assignment
        void operator=(const polarPatchData&);


    // Helper functions

        // return the lats spacing for regular spacing
        scalarList calcLats(const label nlat) const;
        
        // return the lat cell centres
        scalarList calcLatCs(const scalarList&) const;
        
        // calculate the label for the latitude closest to the equator
        label calcJEq(const scalarList& lats) const; 

        // return the lons for regular spacing
        scalarList calcLons(const label nlon) const;
        
        // set the longitudes with unrefinement
        List<scalarList> calcLons
        (
            const scalarList& lats,
            const scalarList& lons,
            const scalarList& urefLats,
            const scalar ratio,
            const label jEq
        ) const;

        // return the lon cell centres
        List<scalarList> calcLonCs(const List<scalarList>& lons) const;


public:

    // Constructors

        //- Construct  for a regular latlon grid and from maxDxLonRsatio
        //  nlon being the maximum number of longitudes for any latitude
        polarPatchData
        (
            const label nlat, 
            const label nlon,
            const scalarList& urefLats,
            const scalar maxDxLonRatio,
            const vector& meshRotation,
            const bool polarCell = true
        );
        
        //- Construct from dictionary
        polarPatchData
        (
            const IOdictionary& earthProperties
        );

    // Access Member Functions

        scalar maxDxLonRatio() const { return maxDxLonRatio_; }

        const scalarList& lats() const { return lats_; }

        const scalarList& latCs() const { return latCs_; }

        const List<scalarList>& lons() const { return lons_; }
        
        const List<scalarList>& lonCs() const { return lonCs_; }
        
        label jEq() const { return jEq_; }

        const scalarList& finestLons() const { return lons_[jEq_]; }

        const scalarList& finestLonCs() const { return lonCs_[jEq_]; }
        
        const vector& meshRotation() const { return meshRotation_; }
        
        const bool& polarCell() const { return polarCell_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
