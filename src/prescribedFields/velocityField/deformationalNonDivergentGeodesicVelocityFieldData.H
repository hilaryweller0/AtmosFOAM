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
    Foam::deformationalNonDivergentGeodesicVelocityFieldData

Description
    Additional data for calculating deformationalNonDivergentGeodesicVelocityField
SourceFiles
    deformationalNonDivergentGeodesicVelocityFieldData.C

\*---------------------------------------------------------------------------*/

#ifndef deformationalNonDivergentGeodesicVelocityFieldData_H
#define deformationalNonDivergentGeodesicVelocityFieldData_H

#include "MeshObject.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class deformationalNonDivergentGeodesicVelocityFieldData Declaration
\*---------------------------------------------------------------------------*/

class deformationalNonDivergentGeodesicVelocityFieldData
:
    public MeshObject<fvMesh, MoveableMeshObject, deformationalNonDivergentGeodesicVelocityFieldData>
{
    // Private Data

        scalar radius_, deformScale_, endTime_;

        vectorField phat_sqrrByT_;
        scalarField Dcos2Lat_;
        scalarField TwoPiSinLat_;

    void calcData();

public:


    // Declare name of the class and its debug switch
    TypeName("deformationalNonDivergentGeodesicVelocityFieldData");

    // Constructors

        //- Construct given an fvMesh and parameters
        explicit deformationalNonDivergentGeodesicVelocityFieldData
        (
            const fvMesh&,
            const scalar radius,
            const scalar deformScale,
            const scalar endTime
        );


    //- Destructor
    ~deformationalNonDivergentGeodesicVelocityFieldData(){}


    // Member Functions

        //- Return references to the data

        const vectorField& phat_sqrrByT() const
        {
            return phat_sqrrByT_;
        }

        const scalarField& Dcos2Lat() const { return Dcos2Lat_; }

        const scalarField& TwoPiSinLat() const { return TwoPiSinLat_; }


        //- Delete the spherical data when the mesh moves
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
