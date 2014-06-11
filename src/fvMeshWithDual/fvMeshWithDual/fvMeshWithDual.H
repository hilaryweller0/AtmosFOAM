/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open So
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::fvMeshWithDual

Description
    Re-reads the face and cell centres, areas and volumes or calculates them
    using spherical geometry and also reads addressing from each
    cell and face onto a patch

SourceFiles
    fvMeshWithDual.C
    fvMeshWithDualGeometry.C
    fvMeshWithDualAddressing.C
    fvMeshWithDualTemplates.C
    fvMeshWithDualSphericalVolGeom.C
    fvMeshWithDualSphericalDistGeom.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshWithDual_H
#define fvMeshWithDual_H

#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "fvBoundaryMesh.H"
#include "surfaceInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvMeshWithDual Declaration
\*---------------------------------------------------------------------------*/

class fvMeshWithDual
:
    public fvMesh
{
    // Private data
    
        // The earthProperties dictionary which holds Omega and earthRadius
        IOdictionary earthProperties_;
    
        //- Four types of geometry
        //-   - CARTESIAN: Cartesian
        //-   - CARTESIANDIST : face and cell centres moved out onto the sphere
        //-                    but distances, areas and vols are Cartesian
        //-   - SPHERICALDIST: face and cell centres moved out onto the sphere,
        //-                    distances defined as great circle distances,
        //-                    areas and vols defined based on distances
        //-   - SPHERICALVOL: everything based on spherical geometry
    public:
        enum geomType{CARTESIAN, CARTESIANDIST, SPHERICALDIST, SPHERICALVOL};
    private:        
        geomType gType_;
        
        //- Read or calculate the cell and face centres
        const Switch readCentres_;
        Switch readCellCentres_;
        Switch readFaceCentres_;

        //- name of the patch for the surface of the sphere
        const word patchName_;
        
        //- patch for the surface of the sphere
        const polyPatch& bottomPatch_;
        
        // Rotation of the sphere
        const dimensionedVector Omega_;
        
        // Magnitude of gravity 
        const dimensionedScalar magg_;
    
        // The axis of rotation of the sphere
        const vector OmegaHat_;
        
        // The reference height
        const dimensionedScalar earthRadius_;
        
        //- Number of levels (or layers) of mesh
        const label nLevels_;

        //- Map from each vertical face onto a patch edge
        mutable IOList<label>* faceToPatchEdgePtr_;

        //- Map from each horizontal face onto a patch face (0 if none
        //- + turning index to patch face addressing)
        mutable IOList<label>* faceToPatchFacePtr_;

        //- Map from each point onto a patch point
        mutable IOList<label>* pointToPatchPointPtr_;
        
        //- Map from each cell onto a patch face
        mutable IOList<label>* cellToPatchFacePtr_;

        //- Reference to the dual mesh (ie the dual of this mesh so it could
        // be the primal of the simulation mesh)
        const fvMeshWithDual* dualMesh_;

        //- Map from each (vertical) face to its dual
        mutable IOList<label>* dualFaceMapPtr_;
        
        //- Orientation of face map from primal to dual
        mutable List<bool>* signFaceMapPtr_;
        
        //- volumes associated with each face corresponding to the owner
        //  and neighbour cell
        mutable surfaceScalarField* faceVolOwnPtr_;
        mutable surfaceScalarField* faceVolNeiPtr_;
        mutable surfaceScalarField* faceVolPtr_;
        
        //- height (depth) of each vertical face (in the vertical)
        mutable surfaceScalarField* depthfPtr_;
        
        //- height (depth) of each cell (in the vertical)
        mutable volScalarField* depthPtr_;
        
        //- area fractions associated with the bottom patch
        mutable List<List<FixedList<scalar,2> > >* fpeAreaFracPtr_;
        mutable scalarListList* facePointAreaFracPtr_;
        
        //- list of points of each face of the patch in same order as above
        mutable labelListList* facePointsPtr_;
        
        //- areas associated with each edge of the patch on each side
        mutable scalarList* edgeAreaOwnPtr_;
        mutable scalarList* edgeAreaNeiPtr_;
        
        //- area of each face of the unit patch
        mutable scalarList* patchFaceAreaPtr_;
        
        //- The diagonal of the H matrix for converting fluxes between primal
        //  and dual on non-orthonogal meshes
        mutable surfaceScalarField* HdiagPtr_;
        
        // The intersection points between primal and dual (near face centres)
        mutable surfaceVectorField* intersectionsPtr_;

        // The volumes associated with each vertical edge of each cell
        mutable scalarListList* cellEdgeVolsPtr_;
        
        // The dual cells associated with each vertical edge of each primal cell
        mutable labelListList* cellEdgeCellsPtr_;

        // Unit vectors in the lon, lat and height directions
        mutable volVectorField* rHatPtr_;
        mutable volVectorField* lonHatPtr_;
        mutable volVectorField* latHatPtr_;
        mutable surfaceVectorField* rHatfPtr_;
        mutable surfaceVectorField* lonHatfPtr_;
        mutable surfaceVectorField* latHatfPtr_;

        // The latitude, longitude and height above reference for each
        // cell and face of the mesh
        mutable volScalarField* lonPtr_;
        mutable volScalarField* latPtr_;
        mutable volScalarField* heightPtr_;
        mutable surfaceScalarField* lonfPtr_;
        mutable surfaceScalarField* latfPtr_;
        mutable surfaceScalarField* heightfPtr_;
        
        // i, j, k dir for each vertical face and ddir (cell centre to centre)
        mutable surfaceVectorField* idirPtr_;
        mutable surfaceVectorField* jdirPtr_;
        mutable surfaceVectorField* kdirPtr_;
        mutable surfaceVectorField* ddirPtr_;
        
    // Private Member Functions

        //- Disallow construct as copy
        fvMeshWithDual(const fvMeshWithDual&);

        //- Disallow assignment
        void operator=(const fvMeshWithDual&);

        //- Read in new face and cell centres if required and calculate
        //-  geometry based on gType_
        void calcGeometry();

        // Read demand driven data
            void readFaceToPatchEdge() const;
            void readFaceToPatchFace() const;
            void readPointToPatchPoint() const;
            void readCellToPatchFace() const;
            void readDualFaceMap() const;
        
        // Make geometric data
            void makesignFaceMap() const;
            void makeSphFaceVols() const;
            void makeFaceDepth() const;
            void makeCellDepth() const;
            void makeSphAreaFracs() const;
            void makeHdiag() const;
            void makeIntersections() const;
            void makeCellEdgeInfo() const;
            void makerHat() const;
            void makelonHat() const;
            void makelatHat() const;
            void makerHatf() const;
            void makelonHatf() const;
            void makelatHatf() const;
            void makelon() const;
            void makelat() const;
            void makeHeight() const;
            void makelonf() const;
            void makelatf() const;
            void makeHeightf() const;
            void makeidir() const;
            void makejdir() const;
            void makekdir() const;
            void makeddir() const;

            void calcSphericalVolGeom() const;
            void calcSphericalDistGeom() const;

        // Helper functions
        const polyPatch& findBottomPatch(const word patchName) const;

public:

    // Static data members

    // Declare name of the class and its debug switch
    ClassName("fvMeshWithDual");

    // Constructors

        //- Construct from IOobject
        explicit fvMeshWithDual
        (
            const IOobject& io,
            const bool readCents = true,
            const geomType gType__ = CARTESIAN
        );

        //- Construct from IOobject and reference to the dual
        explicit fvMeshWithDual
        (
            const IOobject& io,
            const fvMeshWithDual& dualMesh,
            const bool readCents = true
        );


    // Member Functions

        // Access

            //- Type of geometry
            geomType gType() const { return gType_; }
            
            //- Read or calculate the cell and face centres
            Switch readCentres() const { return readCentres_;}

            // The earthProperties dictionary
            const IOdictionary& earthProperties() const
                { return earthProperties_;}
        
            //- name of the patch for the surface of the sphere
            const word& patchName() const { return patchName_; }

            //- patch for the surface of the sphere
            const polyPatch& bottomPatch() const {return bottomPatch_; }

            // Rotation of the sphere
            const dimensionedVector& Omega() const { return Omega_;}
            
            // Magnitude of gravity 
            const dimensionedScalar& magg() const { return magg_;}

            // The axis of rotation of the sphere
            const vector& OmegaHat() const { return OmegaHat_;}
        
            // The reference height
            const dimensionedScalar& earthRadius() const {return earthRadius_;}

            // Number of levels (or layers) of mesh
            label nLevels() const { return nLevels_; }

            //- Map from each vertical face onto a patch edge (-1 if none)
            const IOList<label>& faceToPatchEdge() const;

            //- Map from each horizontal face onto a patch face (0 if none
            //- + turning index to patch face addressing)
            const IOList<label>& faceToPatchFace() const;

            //- Map from horizontal face onto the bottom patch (-1 if none)
            label faceToPatchFace(const label faceI) const;
            
            //- Map from each point onto a patch point
            const IOList<label>& pointToPatchPoint() const;

            //- Map from each cell onto a patch face
            IOList<label> cellToPatchFace() const;

            //- The dual of this mesh
            const fvMeshWithDual& dualMesh() const;

            //- Set the dual of this mesh (modifies face areas of dual for
            //  consistency with primal
            void setDual(const fvMeshWithDual& dm);

            //- Map from each (vertical) face to its dual (-1 if horizontal)
            const IOList<label>& dualFaceMap() const;

            //- Orientation of face map from primal to dual
            const List<bool>& signFaceMap() const;

            //- The sign of the orientation from primal to dual for faceI
            int signMap(const label faceI) const;

            //- volumes associated with each face corresponding to the owner
            const surfaceScalarField& faceVolOwn() const;

            //- volumes associated with each face corresponding to neighbour
            const surfaceScalarField& faceVolNei() const;

            //- total volumes associated with each face
            const surfaceScalarField& faceVol() const;

            //- height (depth) of each vertical face (in the vertical)
            const surfaceScalarField& depthf() const;
            
            //- height (depth) of each cell (in the vertical)
            const volScalarField& depth() const;
            
            //- area fraction of each face of the patch associated with each
            //  patch for each edge
            const List<List<FixedList<scalar,2> > >& fpeAreaFrac() const;

            //- area fraction of each face of patch associated with each point
            const scalarListList& facePointAreaFrac() const;

            //- points of faces of patch in the same order as facePointAreaFrac
            const labelListList& facePoints() const;

            //- areas associated with each edge of the patch on owner side
            const scalarList& edgeAreaOwn() const;

            //- areas associated with each edge of the patch on neighbour side
            const scalarList& edgeAreaNei() const;

            //- area of each face of the unit patch
            const scalarList& patchFaceArea() const;

            //- The diagonal of the H matrix for converting fluxes between
            //  primal and dual on non-orthonogal meshes
            const surfaceScalarField& Hdiag() const
                { if (!HdiagPtr_) makeHdiag(); return *HdiagPtr_;}

            //- The intersection points between primal and dual (near faceCentres)
            const surfaceVectorField& intersections() const;

            //- The volumes associated with each vertical edge of each cell
            const scalarListList& cellEdgeVols() const;
            
            // Dual cells associated with each vertical edge of each primal cell
            const labelListList& cellEdgeCells() const;

            // The latitude, longitude and height above reference for each
            // cell and face of the mesh
            const volScalarField& lon() const;

            const volScalarField& lat() const;

            const volScalarField& height() const;

            const surfaceScalarField& lonf() const;

            const surfaceScalarField& latf() const;

            const surfaceScalarField& heightf() const;

            
            // Unit vectors in the lon, lat and height directions
            const volVectorField& lonHat() const;

            const volVectorField& latHat() const;

            const volVectorField& rHat() const;

            const surfaceVectorField& lonHatf() const;

            const surfaceVectorField& latHatf() const;

            const surfaceVectorField& rHatf() const;

            const surfaceVectorField& idir() const;

            const surfaceVectorField& jdir() const;

            const surfaceVectorField& kdir() const;

            const surfaceVectorField& ddir() const;


            //- Return face deltas as surfaceScalarField
            tmp<surfaceScalarField> magDelta() const;

        // Map fields from one mesh to another
        
            //- Map surfaceField from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;

            //- Map signed surfaceField from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > signedDualMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > signedDualMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;

            //- Map surfaceField of fluxes from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualFluxMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualFluxMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMeshWithDualTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
