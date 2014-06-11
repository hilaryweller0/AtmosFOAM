/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
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

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshWithDual::calcGeometry()
{
    // Based on gType_ and readCentres_ we must calculate and overwrite:
    // faceCentres, cellCentres. faceAreas, cellVolumes, deltaCoeffs

    // First read in face and cell centres if required (and if available)
    if (readCentres_)
    {
        IOobject fCtrsHeader
        (
            "faceCentres", time().findInstance(meshDir(), "points"),
            meshSubDir, *this, IOobject::MUST_READ
        );
        readFaceCentres_ = fCtrsHeader.headerOk();
        if (readFaceCentres_)
        {
            pointIOField fCtrs(fCtrsHeader);
            overrideFaceCentres(fCtrs);
        }
        
        IOobject cCtrsHeader
        (
            "cellCentres", time().findInstance(meshDir(), "points"),
            meshSubDir, *this, IOobject::MUST_READ
        );
        readCellCentres_ = cCtrsHeader.headerOk();
        if (readCellCentres_)
        {
            pointIOField cCtrs(cCtrsHeader);
            overrideCellCentres(cCtrs);
        }
    }
    
    // For CARTESIANDIST face and cell centres are moved out onto the sphere
    // but areas, distances and volumes remain cartesian
    if (gType_ == CARTESIANDIST && !readCentres_)
    {
        FatalErrorIn("fvMeshWithDual::calcGeometry")
            << " gType_ == CARTESIANDIST not yet implemented"
            << exit(FatalError);
    }
    // For SPHERICALVOL: everything based on spherical geometry so no 
    // consistency between primal and dual
    else if (gType_ == SPHERICALVOL)
    {
        calcSphericalVolGeom();
    }
    
    // For SPHERICALDIST, face and cell centres moved out onto the sphere,
    // distances defined as great circle distances, areas and vols defined
    // based on distances. (Requires consistency between primal and dual)

//    else if (gType_ == SPHERICALDIST)
//    {
//        calcSphericalDistGeom();
//    }
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

const Foam::surfaceScalarField& Foam::fvMeshWithDual::faceVolOwn() const
{
    if(!faceVolOwnPtr_) makeSphFaceVols();
    return *faceVolOwnPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshWithDual::faceVolNei() const
{
    if(!faceVolNeiPtr_) makeSphFaceVols();
    return *faceVolNeiPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshWithDual::faceVol() const
{
    if (!faceVolPtr_) makeSphFaceVols();
    return *faceVolPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshWithDual::depthf() const
{
    if (!depthfPtr_) makeFaceDepth();
    return *depthfPtr_;
}


const Foam::volScalarField& Foam::fvMeshWithDual::depth() const
{
    if (!depthPtr_) makeCellDepth();
    return *depthPtr_;
}


const Foam::List<Foam::List<Foam::FixedList<Foam::scalar,2> > >& 
    Foam::fvMeshWithDual::fpeAreaFrac() const
{
    if(!fpeAreaFracPtr_)makeSphAreaFracs();
    return *fpeAreaFracPtr_;
}


const Foam::scalarListList& Foam::fvMeshWithDual::facePointAreaFrac() const
{
    if(!facePointAreaFracPtr_)makeSphAreaFracs();
    return *facePointAreaFracPtr_;
}


const Foam::labelListList& Foam::fvMeshWithDual::facePoints() const
{
    if(!facePointsPtr_) makeSphAreaFracs();
    return *facePointsPtr_;
}


const Foam::scalarList& Foam::fvMeshWithDual::edgeAreaOwn() const
{
    if(!edgeAreaOwnPtr_) makeSphAreaFracs();
    return *edgeAreaOwnPtr_;
}


const Foam::scalarList& Foam::fvMeshWithDual::edgeAreaNei() const
{
    if(!edgeAreaNeiPtr_) makeSphAreaFracs();
    return *edgeAreaNeiPtr_;
}


const Foam::scalarList& Foam::fvMeshWithDual::patchFaceArea() const
{
    if (!patchFaceAreaPtr_) makeSphAreaFracs();
    return *patchFaceAreaPtr_;
}


const Foam::surfaceVectorField& Foam::fvMeshWithDual::intersections() const
{
    if(!intersectionsPtr_) makeIntersections();
    return *intersectionsPtr_;
}


const Foam::scalarListList& Foam::fvMeshWithDual::cellEdgeVols() const
{
    if (!cellEdgeVolsPtr_) makeCellEdgeInfo();
    return *cellEdgeVolsPtr_;
}

const Foam::labelListList& Foam::fvMeshWithDual::cellEdgeCells() const
{
    if (!cellEdgeCellsPtr_) makeCellEdgeInfo();
    return *cellEdgeCellsPtr_;
}


const Foam::volScalarField& Foam::fvMeshWithDual::lon() const
    {if (!lonPtr_) makelon(); return *lonPtr_;}
const Foam::volScalarField& Foam::fvMeshWithDual::lat() const
    { if (!latPtr_) makelat(); return *latPtr_;}
const Foam::volScalarField& Foam::fvMeshWithDual::height() const
    { if (!heightPtr_) makeHeight(); return *heightPtr_;}
const Foam::surfaceScalarField& Foam::fvMeshWithDual::lonf() const
    { if (!lonfPtr_) makelonf(); return *lonfPtr_;}
const Foam::surfaceScalarField& Foam::fvMeshWithDual::latf() const
    { if (!latfPtr_) makelatf(); return *latfPtr_;}
const Foam::surfaceScalarField& Foam::fvMeshWithDual::heightf() const
    { if (!heightfPtr_) makeHeightf(); return *heightfPtr_;}

const Foam::volVectorField& Foam::fvMeshWithDual::lonHat() const
    { if (!lonHatPtr_) makelonHat(); return *lonHatPtr_; }
const Foam::volVectorField& Foam::fvMeshWithDual::latHat() const
    { if (!latHatPtr_) makelatHat(); return *latHatPtr_; }
const Foam::volVectorField& Foam::fvMeshWithDual::rHat() const
    { if (!rHatPtr_) makerHat(); return *rHatPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::lonHatf() const
    { if (!lonHatfPtr_) makelonHatf(); return *lonHatfPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::latHatf() const
    { if (!latHatfPtr_) makelatHatf(); return *latHatfPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::rHatf() const
    { if (!rHatfPtr_) makerHatf(); return *rHatfPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::idir() const
    { if (!idirPtr_) makeidir(); return *idirPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::jdir() const
    { if (!jdirPtr_) makejdir(); return *jdirPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::kdir() const
    { if (!kdirPtr_) makekdir(); return *kdirPtr_; }
const Foam::surfaceVectorField& Foam::fvMeshWithDual::ddir() const
    { if (!ddirPtr_) makeddir(); return *ddirPtr_; }



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMeshWithDual::makeSphFaceVols() const
{
    if (faceVolOwnPtr_ || faceVolNeiPtr_ || faceVolPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makefaceVols")
            << "faceVol already exists" << abort(FatalError);
    }

    faceVolOwnPtr_ = new surfaceScalarField
    (
        IOobject("faceVolOwn", pointsInstance(), *this),
        *this,
        dimensionedScalar("V", dimVol, scalar(0))
    );
    faceVolNeiPtr_ = new surfaceScalarField
    (
        IOobject("faceVolNei", pointsInstance(), *this),
        *this,
        dimensionedScalar("V", dimVol, scalar(0))
    );

    surfaceScalarField& fvo = *faceVolOwnPtr_;
    surfaceScalarField& fvn = *faceVolNeiPtr_;

    // for all vertial faces
    forAll(fvo, faceI)
    {
        const label ie = faceToPatchEdge()[faceI];
        if (ie != -1)
        {
            const label own = owner()[faceI];
            const label nei = neighbour()[faceI];
        
            fvo[faceI] = edgeAreaOwn()[ie]*depth()[own]/3.*
                (3*magSqr(Cf()[faceI]) + 0.25*sqr(depth()[own]));
            fvn[faceI] = edgeAreaNei()[ie]*depth()[nei]/3.*
                (3*magSqr(Cf()[faceI]) + 0.25*sqr(depth()[nei]));
        }
    }

//    // Scale the face vols so that they sum to each cell volume
//    scalarField Vnew(nCells(), scalar(0));
//    forAll(fvo, faceI)
//    {
//        Vnew[owner()[faceI]] += fvo[faceI];
//        Vnew[neighbour()[faceI]] += fvn[faceI];
//    }
//    forAll(fvo.boundaryField(), patchI)
//    {
//        const labelUList& pFaceCells = boundary()[patchI].faceCells();
//        forAll(fvo.boundaryField()[patchI], fi)
//        {
//            Vnew[pFaceCells[fi]] += fvo.boundaryField()[patchI][fi];
//        }
//    }
//    
//    Info << "Face vol difference = " << (Vnew - V())/V() << endl;
    
//    forAll(fvo, faceI)
//    {
//        fvo[faceI] *= V()[owner()[faceI]]/Vnew[owner()[faceI]];
//        fvn[faceI] *= V()[neighbour()[faceI]]/Vnew[neighbour()[faceI]];
//    }
//    forAll(fvo.boundaryField(), patchI)
//    {
//        const labelUList& pFaceCells = boundary()[patchI].faceCells();
//        forAll(fvo.boundaryField()[patchI], fi)
//        {
//            fvo.boundaryField()[patchI][fi] *= V()[pFaceCells[fi]]
//                                              /Vnew[pFaceCells[fi]];
//        }
//    }

    faceVolPtr_ = new surfaceScalarField
    (
        "faceVol",
        faceVolOwn() + faceVolNei()
    );
}

void fvMeshWithDual::makeCellDepth() const
{
    if (depthPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeCellDepth")
            << " depthPtr_ already exists" << exit(FatalError);
    }
    depthPtr_ = new volScalarField
    (
        IOobject
        (
            "depth",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this
        ),
        *this,
        dimensionedScalar("depth", dimLength, scalar(0))
    );
    volScalarField& d = *depthPtr_;
    
    forAll(d, cellI)
    {
        // Find the upper and lower spherical faces of the cell
        label faceUp = -1;
        label faceDown = -1;
        const cell& c = cells()[cellI];
        forAll(c, fi)
        {
            label faceI = c[fi];
            if (faceToPatchFace(faceI) != -1)
            {
                if (faceUp == -1) faceUp = faceI;
                else faceDown = faceI;
            }
        }
        d[cellI] = mag(faceCentres()[faceUp] - faceCentres()[faceDown]);
    }
}

void fvMeshWithDual::makeFaceDepth() const
{
    if (depthfPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeFaceDepth")
            << " depthfPtr_ already exists" << exit(FatalError);
    }
    depthfPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "depthf",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this
        ),
        *this,
        dimensionedScalar("depthf", dimLength, scalar(0))
    );
    surfaceScalarField& d = *depthfPtr_;
    
    forAll(d, faceI)
    {
        // Find the edge of the patch for this (vertical) face
        const label ie = faceToPatchEdge()[faceI];
        if (ie != -1)
        {
            // points at either end of the edge
            const point p0 = unitVector
            (
                bottomPatch_.localPoints()[bottomPatch_.edges()[ie][0]]
            );
            const point p1 = unitVector
            (
                bottomPatch_.localPoints()[bottomPatch_.edges()[ie][1]]
            );
            
            const scalar faceLength = mag(Cf()[faceI])*sphDist(p0,p1);
            d[faceI] = magSf()[faceI]/faceLength;
        }
    }
}

void fvMeshWithDual::makeSphAreaFracs() const
{
    if (fpeAreaFracPtr_ || facePointAreaFracPtr_ || facePointsPtr_
        || edgeAreaOwnPtr_ || edgeAreaNeiPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeSphAreaFracs")
            << "fpeAreaFracPtr_ || facePointAreaFracPtr_ || facePointsPtr_"
            << " already exist" << exit(FatalError);
    }
    
    fpeAreaFracPtr_= new List<List<FixedList<scalar,2> > >(bottomPatch_.size());
    facePointAreaFracPtr_ = new scalarListList(bottomPatch_.size());
    facePointsPtr_ = new labelListList(bottomPatch_.size());
    edgeAreaOwnPtr_ = new scalarList(bottomPatch_.nEdges());
    edgeAreaNeiPtr_ = new scalarList(bottomPatch_.nEdges());
    patchFaceAreaPtr_ = new scalarList(bottomPatch_.size());
    
    List<List<FixedList<scalar,2> > >& fpeAfrac = *fpeAreaFracPtr_;
    scalarListList& fpAfrac = *facePointAreaFracPtr_;
    labelListList& fPts = *facePointsPtr_;
    scalarList& edgeAo = *edgeAreaOwnPtr_;
    scalarList& edgeAn = *edgeAreaNeiPtr_;
    scalarList& fArea  = *patchFaceAreaPtr_;

    // Set the patch edgeCentres from the vertical face centres
    pointField Ce(bottomPatch_.nEdges());
    forAll(Cf(), faceI)
    {
        label ie = faceToPatchEdge()[faceI];
        if (ie != -1) Ce[ie] = intersections()[faceI];
    }
    Ce = unitVector(Ce);

    // Circulate around each point of each face and calculate the proportion
    // of the face area associated with each point and each edge
    forAll(bottomPatch_, faceI)
    {
        // A list of edges for this face (for consistency with the stencil)
        const labelList& faceEdges = bottomPatch_.faceEdges()[faceI];
        
        const face& f = bottomPatch_[faceI];
        const point Ci = unitVector(bottomPatch_.faceCentres()[faceI]);
        
        // The dual point in this face
        //!! Warning, this won't work if either mesh is renumbered!!?!?
        //!! Warning, this won't work if either mesh is renumbered!!?!?
        //const point Ci = unitVector(dualMesh().bottomPatch().points()[faceI]);
        //!! Warning, this won't work if either mesh is renumbered!!?!?
        //!! Warning, this won't work if either mesh is renumbered!!?!?
        
        // Initialise the area fraction arrays and the facePoints array
        fpeAfrac[faceI].setSize(f.size());
        fpAfrac[faceI].setSize(f.size());
        fpAfrac[faceI] = scalar(0);
        fPts[faceI].setSize(f.size());
        
        // total for scaling to make into fractions
        scalar Atot = 0;
        
        // Circulate around face and calculate the half areas for each point
        label prevIp = bottomPatch_.edges()[faceEdges[0]][0];
        if
        (
            prevIp == bottomPatch_.edges()[faceEdges[1]][0]
         || prevIp == bottomPatch_.edges()[faceEdges[1]][1]
        )
        {
            prevIp = bottomPatch_.edges()[faceEdges[0]][1];
        }
        label prevIe = faceEdges.size()-1;
        forAll(faceEdges, iee)
        {
            const label ie = faceEdges[iee];
            label ip = bottomPatch_.edges()[ie][1];
            if (ip == prevIp) ip = bottomPatch_.edges()[ie][0];

            // current and previous points
            const point thisPoint = unitVector(bottomPatch_.localPoints()[ip]);
            const point prevPoint = unitVector(bottomPatch_.localPoints()[prevIp]);
            
            scalar Ap = 0, At = 0;
            if (gType_ == SPHERICALVOL)
            {
                Ap = sphTriSolidAngle(Ce[ie], Ci, prevPoint);
                At = sphTriSolidAngle(Ce[ie], Ci, thisPoint);
                //A = sphTriSolidAngle(Ci, prevPoint, thisPoint);
            }
            else if(gType_ == SPHERICALDIST)
            {
                Ap = sphTriDistAngle(Ce[ie], Ci, prevPoint);
                At = sphTriDistAngle(Ce[ie], Ci, thisPoint);
            }
            else
            {
                FatalErrorIn("makeSphAreaFracs")
                    << "geometry must be one of SPHERICALDIST or SPHERICALVOL "
                    << "not " << gType_ << abort(FatalError);
            }

            fpeAfrac[faceI][iee][0]    = At;
            fpeAfrac[faceI][prevIe][1] = Ap;
            fpAfrac[faceI][iee] += At;
            fpAfrac[faceI][prevIe] += Ap;
            fPts[faceI][iee] = ip;
            Atot += At + Ap;

            if (faceI == bottomPatch_.edgeFaces()[ie][0])
            {
                edgeAo[ie] = At+Ap;
            }
            else
            {
                edgeAn[ie] = At+Ap;
            }
            
            prevIp = ip;
            prevIe = iee;
        }

        // scale to make into fractions
        forAll(faceEdges, iee)
        {
            fpeAfrac[faceI][iee][0] /= Atot;
            fpeAfrac[faceI][iee][1] /= Atot;
            fpAfrac[faceI][iee] /= Atot;
        }
        
        fArea[faceI] = Atot;
    }
}


void fvMeshWithDual::makeHdiag() const
{
    if (HdiagPtr_) FatalErrorIn("fvMeshWithDual::makeHdiag")
                    << "Hdiag already exists" << abort(FatalError);

    HdiagPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "Hdiag", time().findInstance(meshDir(), "points"),
            meshSubDir, *this
        ),
        //dualMesh().signedDualMap(dualMesh().jdir()) & idir()
        ddir() & idir()
    );
    Hdiag().write();
}


void fvMeshWithDual::makeIntersections() const
{
    if (intersectionsPtr_) FatalErrorIn("fvMeshWithDual::makeIntersections")
                        << "intersections already exists" << abort(FatalError);

    intersectionsPtr_ = new surfaceVectorField("intersections", Cf());

    surfaceVectorField& ints = *intersectionsPtr_;

    forAll(ints, faceI)
    {
        label faced = dualFaceMap()[faceI];
        if (faced != -1)
        {
            // Plane of faceI
            const point& p0 = points()[faces()[faceI][0]];
            const point& p1 = points()[faces()[faceI][1]];
            const point& p2 = points()[faces()[faceI][2]];
            plane pFace(p0, p1, p2);
            
            // Plane of faced
            plane dFace
            (
                dualMesh().points()[dualMesh().faces()[faced][0]],
                dualMesh().points()[dualMesh().faces()[faced][1]],
                dualMesh().points()[dualMesh().faces()[faced][2]]
            );
            
            // intersection point
            plane::ray r = pFace.planeIntersect(dFace);
            ints[faceI] = unitVector(r.dir())*sign(r.dir() & Cf()[faceI])
                          *mag(ints[faceI]);
            
            // Warn if intersection does not lie on the primal face
            if
            (
                sign( (ints[faceI] - p0) & (ints[faceI] - p1)) > 0
             && sign( (ints[faceI] - p0) & (ints[faceI] - p2)) > 0
            )
            {
                WarningIn("fvMeshWithDual::makeIntersections")
                    << "primal face " << faceI
                    << " does not overlap with dual face " << faced << endl;
            }
        }
    }
}

void fvMeshWithDual::makeCellEdgeInfo() const
{
    if (cellEdgeVolsPtr_ || cellEdgeCellsPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
             << "cellEdgeVolsPtr or cellEdgeCellsPtr already exists"
             << abort(FatalError);
    }
    
    cellEdgeVolsPtr_ = new scalarListList(nCells());
    cellEdgeCellsPtr_ = new labelListList(nCells());
    
    scalarListList& cellEdgeVols = *cellEdgeVolsPtr_;
    labelListList&  cellEdgeCells = *cellEdgeCellsPtr_;

    // Loop around each primal cell to calculate edgeVols and dual edgeCells
    forAll(cellEdgeVols, cellI)
    {
        scalarList& edgeVols = cellEdgeVols[cellI];
        labelList& edgeCells = cellEdgeCells[cellI];
        
        const cell& c = cells()[cellI];
        // set the number of vertices of the primal cell 2d shape (ie the
        // number of vertical edges)
        const label nVerts = c.size()-2;
        edgeVols.setSize(nVerts);
        edgeCells.setSize(nVerts);
        
        // hash set to store the dual edge cells as we go along
        labelHashSet edgeCellsTmp;
        
        // Also find the radii of the faces above and below cellI
        scalar Rup = -1;
        scalar Rdown = -1;
        
        // Loop around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label faceI = c[fi];
            label faced = dualFaceMap()[faceI];
            if (faced != -1)
            {
                label own = dualMesh().owner()[faced];
                if (!edgeCellsTmp.found(own)) edgeCellsTmp.insert(own);
                label nei = dualMesh().neighbour()[faced];
                if (!edgeCellsTmp.found(nei)) edgeCellsTmp.insert(nei);
            }
            else if (Rup < 0) Rup = mag(Cf()[faceI]);
            else Rdown = mag(Cf()[faceI]);
        }
        // Radii difference for calculating volume
        const scalar Rdiff = Rup > Rdown ? (pow(Rup,3)-pow(Rdown,3))/3.:
                                           (pow(Rdown,3)-pow(Rup,3))/3.;
        
        // Check that edgeCellsTmp is the correct size
        if (edgeCellsTmp.size() != edgeCells.size())
        {
            FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                << "edgeCellsTmp.size() = " << edgeCellsTmp.size()
                << " edgeCells.size() = " << edgeCells.size()
                << " should be the same size" << abort(FatalError);
        }
        edgeCells = edgeCellsTmp.toc();
        
        // Now loop through the edge cells to set the edgeVols
        for(label iv = 0; iv < edgeCells.size(); iv++)
        {
            label cellv = edgeCells[iv];
            // Find the two faces that cellI and cellv have in common
            // (faceI and faceJ) (correspoding to faceD and faceE on dual)
            label faceI = -1;
            label faceJ = -1;
            const cell& cd = dualMesh().cells()[cellv];
            for(label fd = 0; fd < cd.size(); fd++)
            {
                label faced = cd[fd];
                for(label fi = 0; fi < c.size(); fi++)
                {
                    if (faced == dualFaceMap()[c[fi]])
                    {
                        if (faceI == -1)
                        {
                            faceI = c[fi];
                        }
                        else if(faceJ == -1)
                        {
                            faceJ = c[fi];
                        }
                        else
                        {
                            FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                   << "primal and dual cells should not share more than 2 faces"
                            << abort(FatalError);
                        }
                    }
                }
            }
            // Check that faceI and faceJ are both set
            if (faceI == -1 || faceJ == -1)
            {
                FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                    << "2 faces in common not found for primal cell " << cellI
                    << " and dual cell " << cellv << abort(FatalError);
            }
            
            label faceD = dualFaceMap()[faceI];
            label faceE = dualFaceMap()[faceJ];

            // Find one of the points that faceI and faceJ have in common
            label ip = -1;
            const face& facePointsI = faces()[faceI];
            const face& facePointsJ = faces()[faceJ];
            for (label ipi = 0; ipi < facePointsI.size() && ip == -1; ipi++)
            {
                label fpi = facePointsI[ipi];
                for(label ipj = 0; ipj < facePointsJ.size() && ip == -1; ipj++)
                {
                    if (facePointsJ[ipj] == fpi) ip = fpi;
                }
            }
            
            // Find one of the points that faceD and faceE have in common
            label id = -1;
            const face& facePointsD = dualMesh().faces()[faceD];
            const face& facePointsE = dualMesh().faces()[faceE];
            for (label ipi = 0; ipi < facePointsE.size() && id == -1; ipi++)
            {
                label fpi = facePointsD[ipi];
                for(label ipj = 0; ipj < facePointsE.size() && id == -1; ipj++)
                {
                    if (facePointsE[ipj] == fpi) id = fpi;
                }
            }
            
            // We now have primal cell cellI, dual cell cellv, primal faces
            // faceI and faceJ matching dual faces faceD and faceE
            // From these calculate the edgeVols
            //const point& Cp = C()[cellI];
            //const point& Cd = dualMesh().C()[cellv];
            const point& Cp = dualMesh().points()[id];
            const point& Cd = points()[ip];
            const point& Cf0 = intersections()[faceI];
            const point& Cf1 = intersections()[faceJ];
            if (gType_ == SPHERICALVOL)
            {
                edgeVols[iv] = Rdiff*
                (
                    sphTriSolidAngle(Cp,Cd,Cf0)
                  + sphTriSolidAngle(Cp,Cd,Cf1)
                );
            }
            else if(gType_ == SPHERICALDIST)
            {
                edgeVols[iv] = Rdiff*
                (
                    sphTriDistAngle(Cf0, Cp, Cd)
                  + sphTriDistAngle(Cf1, Cp, Cd)
                );
            }
            else
            {
                FatalErrorIn("makeCellEdgeInfo")
                    << "geometry must be one of SPHERICALDIST or SPHERICALVOL "
                    << "not " << gType_ << abort(FatalError);
            }            
        }
//        // Correct so that edgeVols sum to cell vols
//        scalar sumVol = 0;
//        for(label iv = 0; iv < edgeVols.size(); iv++)
//        {
//            sumVol += edgeVols[iv];
//        }
//        for(label iv = 0; iv < edgeVols.size(); iv++)
//        {
//            edgeVols[iv] *= V()[cellI]/sumVol;
//        }
//        Info << "%age difference between primal vols and sum edge vols = "
//             << 100*(sumVol - V()[cellI])/V()[cellI] << nl;
    }
    
//    // Scale edgeVols so that they sum to dual cell vols
//    scalarList dualVols(dualMesh().nCells(), scalar(0));
//    for(label cellI = 0; cellI < nCells(); cellI++)
//    {
//        for(label iv = 0; iv < cellEdgeVols[cellI].size(); iv++)
//        {
//            dualVols[cellEdgeCells[cellI][iv]] += cellEdgeVols[cellI][iv];
//        }
//    }
//    for(label cellI = 0; cellI < nCells(); cellI++)
//    {
//        for(label iv = 0; iv < cellEdgeVols[cellI].size(); iv++)
//        {
//            label cellv = cellEdgeCells[cellI][iv];
//            cellEdgeVols[cellI][iv] *= dualMesh().V()[cellv]/dualVols[cellv];
//        }
//    }
//    Info << "%age difference between dual vols and dual vol sums = "
//         << (dualVols-dualMesh().V())/dualMesh().V()*100 << endl;
}

void fvMeshWithDual::makerHat() const
{
    if (rHatPtr_) FatalErrorIn("fvMeshWithDual::makerHat")
                      << "rHat already exists" << abort(FatalError);
    rHatPtr_ = new volVectorField("rHat", unitVector(C()));
}

void fvMeshWithDual::makelonHat() const
{
    if (lonHatPtr_) FatalErrorIn("fvMeshWithDual::makelonHat")
                        << "lonHat already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonHatPtr_ = new volVectorField("lonHat", unitVector(k ^ rHat()));
}

void fvMeshWithDual::makelatHat() const
{
    if (latHatPtr_) FatalErrorIn("fvMeshWithDual::makelatHat")
                        << "latHat already exists" << abort(FatalError);
    latHatPtr_ = new volVectorField("latHat", unitVector(rHat() ^ lonHat()));
}

void fvMeshWithDual::makerHatf() const
{
    if (rHatfPtr_) FatalErrorIn("fvMeshWithDual::makerHatf")
                        << "rHatf already exists" << abort(FatalError);
    //rHatfPtr_ = new surfaceVectorField("rHatf", unitVector(intersections()));
    rHatfPtr_ = new surfaceVectorField("rHatf", unitVector(Cf()));
}

void fvMeshWithDual::makelonHatf() const
{
    if (lonHatfPtr_) FatalErrorIn("fvMeshWithDual::makelonHatf")
                        << "lonHatf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonHatfPtr_= new surfaceVectorField("lonHatf",unitVector(k^rHatf()));
}

void fvMeshWithDual::makelatHatf() const
{
    if (latHatfPtr_) FatalErrorIn("fvMeshWithDual::makelatHatf")
                        << "latHatf already exists" << abort(FatalError);
    latHatfPtr_= new surfaceVectorField("latHatf",unitVector(rHatf()^lonHatf()));
}

void fvMeshWithDual::makelon() const
{
    if (lonPtr_) FatalErrorIn("fvMeshWithDual::makelon")
                        << "lon already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonPtr_ = new volScalarField
    (
        "lon",
        atan2
        (
            C() & unitVector(vector(0.,1.,0.) - k.y()*k),
            C() & unitVector(vector(1.,0.,0.) - k.x()*k)
        )
    );
}

void fvMeshWithDual::makelat() const
{
    if (latPtr_) FatalErrorIn("fvMeshWithDual::makelat")
                        << "lat already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    latPtr_ = new volScalarField("lon", asin(rHat() & k));
}

void fvMeshWithDual::makeHeight() const
{
    if (heightPtr_) FatalErrorIn("fvMeshWithDual::makeHeight")
                        << "height already exists" << abort(FatalError);
    heightPtr_ = new volScalarField("height", mag(C()) - earthRadius_);
}

void fvMeshWithDual::makelonf() const
{
    if (lonfPtr_) FatalErrorIn("fvMeshWithDual::makelonf")
                        << "lonf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonfPtr_ = new surfaceScalarField
    (
        "lonf",
        atan2
        (
//            intersections() & unitVector(vector(0.,1.,0.) - k.y()*k),
//            intersections() & unitVector(vector(1.,0.,0.) - k.x()*k)
            Cf() & unitVector(vector(0.,1.,0.) - k.y()*k),
            Cf() & unitVector(vector(1.,0.,0.) - k.x()*k)
        )
    );
}

void fvMeshWithDual::makelatf() const
{
    if (latfPtr_) FatalErrorIn("fvMeshWithDual::makelatf")
                        << "latf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    latfPtr_ = new surfaceScalarField("lonf", asin(rHatf() & k));
}

void fvMeshWithDual::makeHeightf() const
{
    if (heightfPtr_) FatalErrorIn("fvMeshWithDual::makeHeightf")
                        << "heightf already exists" << abort(FatalError);
//    heightfPtr_ = new surfaceScalarField("heightf", mag(intersections()) - earthRadius_);
    heightfPtr_ = new surfaceScalarField("heightf", mag(Cf()) - earthRadius_);
}

void fvMeshWithDual::makeidir() const
{
    if (idirPtr_) FatalErrorIn("fvMeshWithDual::makeidir")
                    << "idir already exists" << abort(FatalError);
    idirPtr_ = new surfaceVectorField("idir", unitVector(Sf()));
}

void fvMeshWithDual::makejdir() const
{
    if (jdirPtr_) FatalErrorIn("fvMeshWithDual::makejdir")
                    << "jdir already exists" << abort(FatalError);
    jdirPtr_ = new surfaceVectorField("jdir", kdir() ^ idir());
    
//    surfaceVectorField& j = *jdirPtr_;
//    
//    forAll(j, faceI)
//    {
//        scalar magj = mag(j[faceI]);
//        if (magj < 0.1) j[faceI] = idir()[faceI];
//        else j[faceI] /= magj;
//    }

//    forAll(j.boundaryField(), patchI)
//    {
//        forAll(j.boundaryField()[patchI], patchFace)
//        {
//            vector& ji = j.boundaryField()[patchI][patchFace];
//            scalar magj = mag(ji);
//            if (magj < 0.1)
//            {
//                ji = idir().boundaryField()[patchI][patchFace];
//            }
//            else ji /= magj;
//        }
//    }
}

void fvMeshWithDual::makekdir() const
{
    if (kdirPtr_) FatalErrorIn("fvMeshWithDual::makekdir")
                    << "kdir already exists" << abort(FatalError);
    if
    (
        gType() == SPHERICALDIST || gType() == SPHERICALVOL
     || gType() == CARTESIANDIST
    )
    {
        //kdirPtr_ = new surfaceVectorField("kdir", unitVector(intersections()));
        kdirPtr_ = new surfaceVectorField("kdir", unitVector(Cf()));

        surfaceVectorField& k = *kdirPtr_;
        
        forAll(k, faceI)
        {
            if (faceToPatchEdge()[faceI] == -1) // horizontal face
            k[faceI] = vector(k[faceI].y(), k[faceI].z(), k[faceI].x());
        }

        forAll(k.boundaryField(), patchI)
        {
            forAll(k.boundaryField()[patchI], patchFace)
            {
                vector& ki = k.boundaryField()[patchI][patchFace];
                if
                (
                    faceToPatchEdge()[boundaryMesh()[patchI].start()+patchFace]
                 == -1
                )
                {
                    ki = vector(ki.y(), ki.z(), ki.x());
                }
            }
        }
    }
    else // not spherical
    {
        kdirPtr_ = new surfaceVectorField
        (
            IOobject
            (
                "kdir", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            *this,
            dimensionedVector("kdir", dimless, OmegaHat())
        );
    }
}


void fvMeshWithDual::makeddir() const
{
    if (ddirPtr_) FatalErrorIn("fvMeshWithDual::makeddir")
                    << "ddir already exists" << abort(FatalError);
    ddirPtr_ = new surfaceVectorField("ddir", unitVector(delta()));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
