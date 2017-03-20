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
    Foam::partitionedAtmosphere

Description
    A list of partitions of the atmosphere for the conditinally averaged
    equations. Each partition contains the same number of fluidSpecies in 
    its atmosphere

SourceFiles
    partitionedAtmosphere.C

\*---------------------------------------------------------------------------*/

#ifndef partitionedAtmosphere_H
#define partitionedAtmosphere_H

#include "partition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class partitionedAtmosphere Declaration
\*---------------------------------------------------------------------------*/

class partitionedAtmosphere
:
    public PtrList<partition>
{
    // Private data
    
        //- Total density over all partitions
        volScalarField rho_;

        //- Rate of change of density over all partitions
        volScalarField dRhodt_;
        
        //- Theta density averaged over all partitions
        volScalarField theta_;
        
        //- Velocity at faces
        surfaceVectorField Uf_;
        
        //- Flux at faces
        surfaceScalarField flux_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        partitionedAtmosphere(const partitionedAtmosphere&);

        //- Disallow default bitwise assignment
        void operator=(const partitionedAtmosphere&);


public:

    // Constructors

        //- Construct
        partitionedAtmosphere
        (
            const wordList& partitionNames,
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict
        );


    //- Destructor
    ~partitionedAtmosphere();


    // Member Functions

    // Access
        //- Total density over all partitions
        volScalarField& rho() { return rho_;}
        const volScalarField& rho() const { return rho_;}

        //- Rate of change of density over all partitions
        volScalarField& dRhodt() { return dRhodt_;}
        const volScalarField& dRhodt() const { return dRhodt_;}
        
        //- Theta density averaged over all partitions
        volScalarField& theta() { return theta_; }
        const volScalarField& theta() const { return theta_; }
        
        //- Velocity at faces
        surfaceVectorField& Uf() { return Uf_;}
        const surfaceVectorField& Uf() const { return Uf_;}
        
        //- Flux at faces
        surfaceScalarField& flux() { return flux_;}
        const surfaceScalarField& flux() const { return flux_;}


    //- Update and return fields summed over partitions
        //- Total density
        volScalarField& sumDensity();
        
//        //- Equivalent potential temperature
//        tmp<volScalarField> thetae() const;
//        
        //- Update theta from partitions
        volScalarField& updateTheta();
        
        //- Update Uf from partitions
        surfaceVectorField& updateUf();

        //- Update flux from partitions
        surfaceScalarField& updateFlux();
        
        //- Update volume fractions in all the partitions
        void updateSigmas(const volScalarField& p);

        //- Density weighted total gas constant
        tmp<volScalarField> rhoR() const;
        
        //- Total volume fraction of liquid
        tmp<volScalarField> volLiquid() const;

        // Enxer from the equation of state
        tmp<volScalarField> ExnerFromState();

        // Enxer from the equation of state
        volScalarField& ExnerFromState(volScalarField& Exner);

        //- The coefficient in from of grad(Exner) for the Helmholtz eqn
        void setGradPcoeff(surfaceScalarField& gradPcoeff);

        // Output
        void write();

        // Input (re-read all fields)
        void readUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
