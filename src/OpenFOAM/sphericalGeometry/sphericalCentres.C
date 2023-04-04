/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "sphericalCentres.H"

void Foam::sphericalCentres(const primitiveMesh& mesh)
{
    Info << "Modifying cell centres and face centres for spherical geometry"
         << endl;

    // Centres to change
    vectorField& Cf = const_cast<vectorField&>(mesh.faceCentres());
    vectorField& C = const_cast<vectorField&>(mesh.cellCentres());
    // Need to change volumes to be consistent
    scalarField& V = const_cast<scalarField&>(mesh.cellVolumes());

    Info << "C[0] = " << C[0] << "\nCf[0] = " << Cf[0] << "\nV[0] = "
         << V[0] << endl;
         
    // Copy of old face centres needed for calculating new cell centres
    const vectorField Cfold = mesh.faceCentres();

    // Other mesh data
    const vectorField& Sf(mesh.faceAreas());
    const pointField& p = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    Info << "Cf[0] - C[own[0]] = " << Cf[0] - C[own[0]] << endl;

    // Loop through all faces and modify face centres
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
            // Add area of the triangle between the edge and the cell centre
            scalar At = mag((nextPoint - Cf[facei])^(p[f[ip]] - Cf[facei]));
            A += At;
            // Add average radius for this triangle
            Ar += 0.5*(mag(nextPoint) + mag(p[f[ip]]))*At;
        }
        
        // Push the face centre onto the sphere
        scalar magCf = mag(Cf[facei]);
        if (magCf > SMALL)
        {
            Cf[facei] *= Ar/A/magCf;
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
            scalar Vp = mag((Cfold[facei] - C[celli]) & Sf[facei]);
            Vc += Vp;
            // Add face radius for this face
            Vr += mag(Cf[facei])*Vp;
        }
        
        // Push the cell centre onto the sphere
        scalar magC = mag(C[celli]);
        if (magC > SMALL)
        {
            C[celli] *= Vr/Vc/magC;
        }
    }
    
    // Modify the cell volumes for consistency with face and cell centres
    V = 0;
    forAll(own, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf[facei] & (Cf[facei] - C[own[facei]]);

        // Accumulate face-pyramid volume
        V[own[facei]] += pyr3Vol;
    }

    forAll(nei, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf[facei] & (C[nei[facei]] - Cf[facei]);

        // Accumulate face-pyramid volume
        V[nei[facei]] += pyr3Vol;
    }

    V *= (1.0/3.0);

    Info << "C[0] = " << C[0] << "\nCf[0] = " << Cf[0] << "\nV[0] = "
         << V[0] << endl;
    Info << "Cf[0] - C[own[0]] = " << Cf[0] - C[own[0]] << endl;
}

// ************************************************************************* //
