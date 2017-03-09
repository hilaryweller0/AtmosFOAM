/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::partition

Description
    A partition of the atmosphere for solving the conditionally averaged
    Euler equations conditioned based on the presence of buoyant convection.
    A partition consists of an atmosphere, a temperature field, a volume
    fraction for the partition (sigma) and a velocity field

SourceFiles
    partition.C

\*---------------------------------------------------------------------------*/

#ifndef partition_H
#define partition_H

#include "atmosphere.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class partition Declaration
\*---------------------------------------------------------------------------*/

class partition
:
    public atmosphere
{
    // Private data

        //- Name of this partition
        const word partitionName_;
        
        //- Volume fraction of this partition
        volScalarField sigma_;
        
        //- Density of the partition
        volScalarField rho_;
        
        //- Temperature of this partition
        volScalarField T_;

        //- Potential temperature of this partition
        volScalarField theta_;
        
        //- Velocity of the partition
        surfaceVectorField Uf_;
        
        //- Flux of the partition
        surfaceScalarField flux_;
        
        //- Rate of change of rhoSigma
        volScalarField drhoSigmadt_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        partition(const partition&);

        //- Disallow default bitwise assignment
        void operator=(const partition&);


public:

    // Constructors

        //- Construct from components, dictionaries and IOobjects
        partition
        (
            const word& partitionName__,
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict
        );

        //- Construct from components, dictionaries and IOobjects and Exner
        partition
        (
            const word& partitionName__,
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict,
            const volScalarField& Exner
        );


    //- Destructor
    ~partition();


    // Member Functions

        // Access
        const word& partitionName() const { return partitionName_; }
        
        //- Volume fraction of this partition
        const volScalarField& sigma() const { return sigma_; }
        volScalarField& sigma() { return sigma_; }
        
        //- Density of the partition
        const volScalarField& rho() const { return rho_; }
        volScalarField& rho() { return rho_; }

        //- Temperature of this partition
        const volScalarField& T() const { return T_; }
        volScalarField& T() { return T_; }
        
        //- Potential temperature of this partition
        const volScalarField& theta() const { return theta_; }
        volScalarField& theta() { return theta_; }
        
        //- Velocity of the partition
        const surfaceVectorField& Uf() const { return Uf_; }
        surfaceVectorField& Uf() { return Uf_; }
        
        //- Flux of the partition
        const surfaceScalarField& flux() const { return flux_; }
        surfaceScalarField& flux() { return flux_; }
        
        //- Rate of change of rhoSigma
        const volScalarField& drhoSigmadt() const { return drhoSigmadt_; }
        volScalarField& drhoSigmadt() { return drhoSigmadt_; }

        //- Update theta and T from Exner
        void updateThetaT(const volScalarField& Exner);

        //- Update rho from the phases
        volScalarField& sumDensity() const
        {
            rho_ = atmosphere::sumDensity();
            return rho_;
        }

        // Write
        void write();

        // Input (re-read densities, liquid fractions, temperature and sigma)
        void readUpdate(const volScalarField& Exner);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
