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
    Foam::fluidSpecie

Description
    Contains a perfect gas and liquid phase and constant and functions
    describing phase changes between them

SourceFiles
    fluidSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSpecie_H
#define fluidSpecie_H

#include "perfectGasPhase.H"
#include "liquidPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluidSpecie Declaration
\*---------------------------------------------------------------------------*/

class fluidSpecie
{
    // Private data

        //- Name of this species
        const word name_;

        //- gas and liquid phases
        perfectGasPhase gas_;
        liquidPhase liquid_;
        
        //- Constants in equations for phase changes
        const dimensionedScalar Lv0_;
        const dimensionedScalar pvs0_;

        //- Condensation of this phase
        volScalarField condensation_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluidSpecie(const fluidSpecie&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSpecie&);


public:

    // Constructors

        //- Construct from Istreams, dictionaries and IOobjects
        fluidSpecie
        (
            const word& name__,
            const IOobject& gasIO,
            const IOobject& liquidIO,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~fluidSpecie();


    // Member Functions

        // Access

            //- gas and liquid phases
            const word& name() { return name_; }
            perfectGasPhase& gas() { return gas_; }
            liquidPhase& liquid() { return liquid_; }
            
            // Constants in equations for phase changes
            const dimensionedScalar& Lv0() { return Lv0_; }
            const dimensionedScalar& pvs0() { return pvs0_; }

            //- Latent heat of vaporisation
            tmp<volScalarField> latentHeat(const volScalarField& T);

            //- Saturation vapour pressure
            tmp<volScalarField> pSat(const volScalarField& T);

            //- Condensation in one time interval
            const volScalarField& condensation() { return condensation_; }
            const volScalarField& updateCondensation(const volScalarField& T);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
