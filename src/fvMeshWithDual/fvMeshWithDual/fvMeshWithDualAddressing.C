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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::polyPatch& Foam::fvMeshWithDual::findBottomPatch
(
    const word patchName
) const
{
    label patchID = boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        return boundaryMesh().last();
    }
    return boundaryMesh()[patchID];
}


void Foam::fvMeshWithDual::readFaceToPatchEdge() const
{
    if (faceToPatchEdgePtr_) FatalErrorIn("fvMeshWithDual::readFaceToPatchEdge")
                         << "faceToPatchEdge already exists" << abort(FatalError);
    faceToPatchEdgePtr_ = new IOList<label>
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readFaceToPatchFace() const
{
    if (faceToPatchFacePtr_) FatalErrorIn("fvMeshWithDual::readFaceToPatchFace")
                         << "faceToPatchFace already exists" << abort(FatalError);
    faceToPatchFacePtr_ = new IOList<label>
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readPointToPatchPoint() const
{
    if(pointToPatchPointPtr_)FatalErrorIn("fvMeshWithDual::readPointToPatchPoint")
                        << "pointToPatchPoint already exists" << abort(FatalError);
    pointToPatchPointPtr_ = new IOList<label>
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readCellToPatchFace() const
{
    if (cellToPatchFacePtr_) FatalErrorIn("fvMeshWithDual::readCellToPatchFace")
                         << "cellToPatchFace already exists" << abort(FatalError);
    cellToPatchFacePtr_ = new IOList<label>
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readDualFaceMap() const
{
    if (dualFaceMapPtr_) FatalErrorIn("fvMeshWithDual::readDualFaceMap")
                         << "dualFaceMap already exists" << abort(FatalError);
    dualFaceMapPtr_ = new IOList<label>
    (
        IOobject
        (
            "dualFaceMap",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}


void Foam::fvMeshWithDual::makesignFaceMap() const
{
    if (signFaceMapPtr_) FatalErrorIn("fvMeshWithDual::makesignFaceMap")
                           << "signFaceMap already exists" << abort(FatalError);
    signFaceMapPtr_ = new List<bool>(dualFaceMap().size());
    List<bool>& sfm = *signFaceMapPtr_;

    forAll(dualFaceMap(), faceI)
    {
        if (dualFaceMap()[faceI] != -1)
        {
            label faced = dualFaceMap()[faceI];
            sfm[faceI] = sign
            (
                (faceAreas()[faceI] ^ dualMesh().faceAreas()[faced])
              & faceCentres()[faceI]
            ) == 1;
        }
    }
}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::IOList<Foam::label>& Foam::fvMeshWithDual::faceToPatchEdge() const
{
    if(!faceToPatchEdgePtr_) readFaceToPatchEdge();
    return *faceToPatchEdgePtr_;
}


const Foam::IOList<Foam::label>& Foam::fvMeshWithDual::faceToPatchFace() const
{
    if (!faceToPatchFacePtr_) readFaceToPatchFace();
    return *faceToPatchFacePtr_;
}


Foam::label Foam::fvMeshWithDual::faceToPatchFace(const label faceI) const
{
    label fi = faceToPatchFace()[faceI];
    if (fi == 0) return -1;
    fi = mag(fi)-1;
    label patchID = boundaryMesh().whichPatch(fi);
    return fi - boundaryMesh()[patchID].start();
}


const Foam::IOList<Foam::label>& Foam::fvMeshWithDual::pointToPatchPoint() const
{
    if (!pointToPatchPointPtr_) readPointToPatchPoint();
    return *pointToPatchPointPtr_;
}


Foam::IOList<Foam::label> Foam::fvMeshWithDual::cellToPatchFace() const
{
    if (!cellToPatchFacePtr_) readCellToPatchFace();
    return *cellToPatchFacePtr_;
}


const Foam::IOList<Foam::label>& Foam::fvMeshWithDual::dualFaceMap() const
{
    if (!dualFaceMapPtr_) readDualFaceMap();
    return *dualFaceMapPtr_;
}


const Foam::List<bool>& Foam::fvMeshWithDual::signFaceMap() const
{
    if(!signFaceMapPtr_) makesignFaceMap();
    return *signFaceMapPtr_;
}


int Foam::fvMeshWithDual::signMap(const label faceI) const
{
    return 2*signFaceMap()[faceI]-1;
}

// ************************************************************************* //
