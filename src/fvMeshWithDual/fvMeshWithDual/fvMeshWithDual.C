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
#include "scalarIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshWithDual, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshWithDual::fvMeshWithDual
(
    const IOobject& io,
    const bool readCents,
    const fvMeshWithDual::geomType gType__
)
:
    fvMesh(io),
    earthProperties_
    (
        IOobject
        (
            "earthProperties", time().constant(), *this, IOobject::READ_IF_PRESENT
        )
    ),
    gType_(gType__),
    readCentres_(readCents),
    readCellCentres_(readCents),
    readFaceCentres_(readCents),
    patchName_(earthProperties_.lookupOrDefault<word>("surfacePatch", "originalPatch")),
    bottomPatch_(findBottomPatch(patchName_)),
    Omega_(earthProperties_.lookupOrDefault<dimensionedVector>("Omega", dimensionedVector("O", dimless, vector(0.,0.,0.)))),
    magg_(earthProperties_.lookupOrDefault<dimensionedScalar>("magg", dimensionedScalar("magg", dimless, scalar(0)))),
    OmegaHat_(unitVector(Omega_.value())),
    earthRadius_(earthProperties_.lookupOrDefault<dimensionedScalar>("earthRadius", dimensionedScalar("earthRadius", dimLength, scalar(1)))),
    nLevels_(nCells()/(max(bottomPatch_.size(),1))),
    faceToPatchEdgePtr_(NULL),
    faceToPatchFacePtr_(NULL),
    pointToPatchPointPtr_(NULL),
    cellToPatchFacePtr_(NULL),
    dualMesh_(NULL),
    dualFaceMapPtr_(NULL),
    signFaceMapPtr_(NULL),
    faceVolOwnPtr_(NULL),
    faceVolNeiPtr_(NULL),
    faceVolPtr_(NULL),
    depthfPtr_(NULL),
    depthPtr_(NULL),
    fpeAreaFracPtr_(NULL),
    facePointAreaFracPtr_(NULL),
    facePointsPtr_(NULL),
    edgeAreaOwnPtr_(NULL),
    edgeAreaNeiPtr_(NULL),
    patchFaceAreaPtr_(NULL),
    HdiagPtr_(NULL),
    intersectionsPtr_(NULL),
    cellEdgeVolsPtr_(NULL),
    cellEdgeCellsPtr_(NULL),
    rHatPtr_(NULL),
    lonHatPtr_(NULL),
    latHatPtr_(NULL),
    rHatfPtr_(NULL),
    lonHatfPtr_(NULL),
    latHatfPtr_(NULL),
    lonPtr_(NULL),
    latPtr_(NULL),
    heightPtr_(NULL),
    lonfPtr_(NULL),
    latfPtr_(NULL),
    heightfPtr_(NULL),
    idirPtr_(NULL),
    jdirPtr_(NULL),
    kdirPtr_(NULL),
    ddirPtr_(NULL)
{
    if
    (
        gType_ != CARTESIAN && gType_ != CARTESIANDIST
     && gType_ != SPHERICALDIST && gType_ != SPHERICALVOL
    )
    {
        FatalErrorIn("fvMeshWithDual::fvMeshWithDual") << "geomType must be one of fvMeshWithDual::CARTESIAN, fvMeshWithDual::CARTESIANDIST, fvMeshWithDual::SPHERICALDIST, fvMeshWithDual::SPHERICALVOL. Not " << gType_
            << abort(FatalError);
    }

    Info << "Creating mesh with geometry " << gType_ << endl;
    calcGeometry();
    polyMesh::boundary_.updateMesh();
    polyMesh::boundary_.calcGeometry();
    
    if (gType_ == CARTESIAN)
    {
        surfaceScalarField dc(deltaCoeffs());
        surfaceInterpolation::overrideDeltaCoeffs(dc);
    }
}

Foam::fvMeshWithDual::fvMeshWithDual
(
    const IOobject& io,
    const fvMeshWithDual& dualMesh__,
    const bool readCents
)
:
    fvMesh(io),
    earthProperties_
    (
        IOobject
        (
            "earthProperties", time().constant(), *this, IOobject::READ_IF_PRESENT
        )
    ),
    gType_(dualMesh__.gType()),
    readCentres_(readCents),
    readCellCentres_(readCents),
    readFaceCentres_(readCents),
    patchName_(earthProperties_.lookupOrDefault<word>("surfacePatch", "originalPatch")),
    bottomPatch_(findBottomPatch(patchName_)),
    Omega_(earthProperties_.lookupOrDefault<dimensionedVector>("Omega", dimensionedVector("O", dimless, vector(0.,0.,0.)))),
    magg_(earthProperties_.lookupOrDefault<dimensionedScalar>("magg", dimensionedScalar("magg", dimless, scalar(0)))),
    OmegaHat_(unitVector(Omega_.value())),
    earthRadius_(earthProperties_.lookupOrDefault<dimensionedScalar>("earthRadius", dimensionedScalar("earthRadius", dimLength, scalar(1)))),
    nLevels_(nCells()/(max(bottomPatch_.size(),1))),
    faceToPatchEdgePtr_(NULL),
    faceToPatchFacePtr_(NULL),
    pointToPatchPointPtr_(NULL),
    cellToPatchFacePtr_(NULL),
    dualMesh_(&dualMesh__),
    dualFaceMapPtr_(NULL),
    signFaceMapPtr_(NULL),
    faceVolOwnPtr_(NULL),
    faceVolNeiPtr_(NULL),
    faceVolPtr_(NULL),
    depthfPtr_(NULL),
    depthPtr_(NULL),
    fpeAreaFracPtr_(NULL),
    facePointAreaFracPtr_(NULL),
    facePointsPtr_(NULL),
    edgeAreaOwnPtr_(NULL),
    edgeAreaNeiPtr_(NULL),
    patchFaceAreaPtr_(NULL),
    HdiagPtr_(NULL),
    intersectionsPtr_(NULL),
    cellEdgeVolsPtr_(NULL),
    cellEdgeCellsPtr_(NULL),
    rHatPtr_(NULL),
    lonHatPtr_(NULL),
    latHatPtr_(NULL),
    rHatfPtr_(NULL),
    lonHatfPtr_(NULL),
    latHatfPtr_(NULL),
    lonPtr_(NULL),
    latPtr_(NULL),
    heightPtr_(NULL),
    lonfPtr_(NULL),
    latfPtr_(NULL),
    heightfPtr_(NULL),
    idirPtr_(NULL),
    jdirPtr_(NULL),
    kdirPtr_(NULL),
    ddirPtr_(NULL)
{
    Info << "Creating mesh and dual mesh with geometry " << gType_ << endl;
    calcGeometry();
    if (gType_ == SPHERICALDIST)
    {
        calcSphericalDistGeom();
    }
    polyMesh::boundary_.updateMesh();
    polyMesh::boundary_.calcGeometry();
}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



const Foam::fvMeshWithDual& Foam::fvMeshWithDual::dualMesh() const
{
    if (!dualMesh_)
    {
        FatalErrorIn("fvMeshWithDual::dualMesh")
            << "dual mesh not set" << abort(FatalError);
    }
    return *dualMesh_;
}


void Foam::fvMeshWithDual::setDual(const fvMeshWithDual& dm)
{
    dualMesh_ = &dm;
    
    if (gType_ == SPHERICALDIST)
    {
        calcSphericalDistGeom();
    }
    else 
    if (gType_ == SPHERICALVOL)
    {
        // re-calculate and override deltaCoeffs    
        surfaceScalarField deltaCoeffs(0.5*magSf()/faceVol());
        surfaceInterpolation::overrideDeltaCoeffs(deltaCoeffs);
    }
}


// ************************************************************************* //
