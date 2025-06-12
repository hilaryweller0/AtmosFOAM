/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "sphericalMeshData.H"
#include "polarPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphericalMeshData, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::sphericalMeshData::sphericalMeshData
(
    const fvMesh& mesh,
    const scalar earthRadius__
)
:
    DemandDrivenMeshObject<fvMesh, Foam::MoveableMeshObject, sphericalMeshData>(mesh),
    Cf_(mesh.faceCentres()),
    C_(mesh.cellCentres()),
    V_(mesh.nCells()),
    Sf_(mesh.faceAreas()),
    earthRadius_(earthRadius__),
    CfLonLatz_(mesh.nFaces()),
    CLonLatz_(mesh.nCells()),
    pointsLatLonz_(mesh.nPoints())
{
    calcSphericalMeshData();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::sphericalMeshData::~sphericalMeshData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphericalMeshData::calcSphericalMeshData()
{
    if (debug)
    {
        InfoInFunction << "Calculating spherical geometry" << endl;
    }

    const fvMesh& mesh = this->mesh();

    // References to Cartesian mesh data
    const vectorField& CfC = mesh.faceCentres();
    const vectorField& Sf(mesh.faceAreas());
    const pointField& p = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Loop through all faces and calculate face centres on the sphere
    forAll(faces, facei)
    {
        const labelList& f = faces[facei];
        label nPoints = f.size();
        
        // Sum the area weighted face radius and the area
        scalar Ar = 0;
        scalar A = 0;
        
        // Loop around the face to sum the area weighted face radius and area
        for(label ip = 0; ip < nPoints; ip++)
        {
            const point& nextPoint = p[f[(ip+1)%nPoints]];
            // Add area of the triangle between the edge and the face centre
            scalar At = mag((nextPoint - CfC[facei])^(p[f[ip]] - CfC[facei]));
            A += At;
            // Add average radius for this triangle
            Ar += 0.5*(mag(nextPoint) + mag(p[f[ip]]))*At;
        }

        // Push the face centre onto the sphere
        scalar magCf = mag(CfC[facei]);
        if (magCf > SMALL)
        {
            Cf_[facei] *= Ar/A/magCf;
        }
    }

    // Loop through all cells and modify cell centres (using old face centres)
    forAll(cells, celli)
    {
        const labelList& c = cells[celli];
        label nFaces = c.size();
        
        // Sum the volume weighted face radius and the volume
        scalar Vr = 0;
        scalar Vc = 0;
        
        // Loop over the faces of the cell to sum the volume weighted radius
        for(label fi = 0; fi < nFaces; fi++)
        {
            label facei = c[fi];
            // The volume of the prism from the old face to the old centre
            scalar Vp = mag((CfC[facei] - C_[celli]) & Sf[facei]);
            Vc += Vp;
            // Add face radius for this face
            Vr += mag(Cf_[facei])*Vp;
        }
        
        // Push the cell centre onto the sphere
        scalar magC = mag(C_[celli]);
        if (magC > SMALL)
        {
            C_[celli] *= Vr/Vc/magC;
        }
    }
    
    // Modify the cell volumes for consistency with face and cell centres
    V_ = 0;
    forAll(own, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf[facei] & (Cf_[facei] - C_[own[facei]]);

        // Accumulate face-pyramid volume
        V_[own[facei]] += pyr3Vol;
    }

    forAll(nei, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf[facei] & (C_[nei[facei]] - Cf_[facei]);

        // Accumulate face-pyramid volume
        V_[nei[facei]] += pyr3Vol;
    }

    V_ *= (1.0/3.0);
    
    // Calculating CfLonLatz_ and  CLonLatz_
    CfLonLatz_.replace(2, mag(Cf_));
    CfLonLatz_.replace(0, atan2(Cf_.component(1), Cf_.component(0)));
    CfLonLatz_.replace
    (
        1,
        asin(max(min(Cf_.component(2)/CfLonLatz_.component(2), 1.), -1.))
    );
    CfLonLatz_.replace(2, CfLonLatz_.component(2) - earthRadius_);
    CLonLatz_.replace(2, mag(C_));
    CLonLatz_.replace(0, atan2(C_.component(1), C_.component(0)));
    CLonLatz_.replace
    (
        1,
        asin(max(min(C_.component(2)/CLonLatz_.component(2), 1.), -1.))
    );
    CLonLatz_.replace(2, CLonLatz_.component(2) - earthRadius_);
    
    // Calculating pointsLatLonz_
    forAll(pointsLatLonz_, ip)
    {
        polarPoint pp = convertToPolar(p[ip], 180, earthRadius_);
        pointsLatLonz_[ip] = point(pp[0], pp[1], pp[2]);
    }
    
    if (debug)
    {
        InfoInFunction
            << "Finished calculating spherical geometry" << endl;
    }
}


bool Foam::sphericalMeshData::movePoints()
{
    calcSphericalMeshData();
    return true;
}


// ************************************************************************* //
