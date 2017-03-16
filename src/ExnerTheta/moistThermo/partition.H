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

#include "baseAtmosphere.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class partition Declaration
\*---------------------------------------------------------------------------*/

class partition
:
    public baseAtmosphere
{
    // Private data

        //- Name of this partition
        const word partitionName_;
        
        //- Volume fraction of the partition
        volScalarField sigma_;

        //- Density of the partition multiplied by the volume ratio
        volScalarField sigmaRho_;
        
        //- Potential temperature of this partition
        volScalarField theta_;
        
        //- Temperature of this partition
        volScalarField T_;

        //- Velocity of the partition
        surfaceVectorField Uf_;
        
        //- Flux of the partition (sigmaRho * Uf).Sf
        surfaceScalarField flux_;

        //- Volume of liquid in the partition
        volScalarField sigmaVolLiquid_;

        //- Rate of change of density of the partition
        volScalarField dSigmaRhodt_;

        //- Rate of change of rhoTheta of the partition
        volScalarField dSigmaRhoThetadt_;
        
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


    //- Destructor
    ~partition();


    // Member Functions

    // Access
        const word& partitionName() const { return partitionName_; }

        //- Volume fraction of the partition
        volScalarField& sigma() { return sigma_;}
        const volScalarField& sigma() const { return sigma_;}

        //- Density of the partition
        const volScalarField& sigmaRho() const { return sigmaRho_; }
        volScalarField& sigmaRho() { return sigmaRho_; }

        //- Rate of change of density of the partition
        const volScalarField& dSigmaRhodt() const { return dSigmaRhodt_; }
        volScalarField& dSigmaRhodt() { return dSigmaRhodt_; }

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
        
        //- Volume of liquid in the partition
        volScalarField& sigmaVolLiquid() { return sigmaVolLiquid_;}
        const volScalarField& sigmaVolLiquid() const { return sigmaVolLiquid_;}

        //- Rate of change of rhoTheta of the partition
        volScalarField& dSigmaRhoThetadt() { return dSigmaRhoThetadt_;}
        const volScalarField& dSigmaRhoThetadt() const
             { return dSigmaRhoThetadt_;}

    // Update functions

        //- Update rho from the phases
        volScalarField& sumDensity();

        //- Update volLiquied from the phases
        volScalarField& sumVolLiquid();
        
        //- Calculate sigma from the total pressure
        volScalarField& updateSigma(const volScalarField& p);

    // Constant member functions

        // Source term in the theta equation
        tmp<volScalarField> thetaSource() const;
        
        // Enxer from the equation of state
        tmp<volScalarField> exnerFromState() const;
        
        // Write
        void write();

        // Input (re-read densities, liquid fractions and theta)
        void readUpdate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
