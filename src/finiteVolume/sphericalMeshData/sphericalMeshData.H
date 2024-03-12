/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::sphericalMeshData

Description
    Additional mesh data for spherical geometry:
    Cf: faceCentres are pushed out onto the same sphere as the verticies
    C:  cellCentres are calculated from faceCentres on the sphere
    V:  volumes are consistent with Cf and C
    Sf: face areas are consistent with Cf and verticies (but not using 
        spherical triangles)
SourceFiles
    sphericalMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalMeshData_H
#define sphericalMeshData_H

#include "DemandDrivenMeshObject.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sphericalMeshData Declaration
\*---------------------------------------------------------------------------*/

class sphericalMeshData
:
    public DemandDrivenMeshObject<fvMesh, MoveableMeshObject, sphericalMeshData>
{
    // Private Data

        //- Alternative mesh data for spherical geometry
        vectorField Cf_;
        vectorField C_;
        scalarField V_;
        vectorField Sf_;
        
        //- Additional mesh data for spherical geomtetry
        scalar earthRadius_;
        vectorField CfLonLatz_;
        vectorField CLonLatz_;
        pointField pointsLatLonz_;

    // Private Member Functions

        //- Construct spherical geometry
        void calcSphericalMeshData();


public:

    // Static data members

        static scalar constexpr degToRad = M_PI/180.;
        static scalar constexpr radToDeg = 180./M_PI;

    // Declare name of the class and its debug switch
    TypeName("sphericalMeshData");


    // Constructors

        //- Construct given an fvMesh and earthRadius
        explicit sphericalMeshData
        (
            const fvMesh&,
            const scalar earthRadius__ = 0
        );


    //- Destructor
    virtual ~sphericalMeshData();


    // Member Functions

        //- Return reference to the face centres
        const vectorField& faceCentres() const
        {
            return Cf_;
        }

        //- Return reference to the cell centres
        const vectorField& cellCentres() const { return C_; }

        //- Return reference to the cell volumes
        const scalarField& V() const { return V_; }

        //- Return reference to the face area vectors
        const vectorField& faceAreas() const { return Sf_; }

        scalar earthRadius() const { return earthRadius_; }
        const vectorField& CfLonLatz() const { return CfLonLatz_; }
        const vectorField& CLonLatz() const { return CLonLatz_; }
        const pointField& pointsLatLonz() const { return pointsLatLonz_; }

        //- Delete the spherical data when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
