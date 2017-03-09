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

        partitionedAtmosphere
        (
            const wordList& partitionNames,
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict,
            const volScalarField& Exner
        );


    //- Destructor
    ~partitionedAtmosphere();


    // Member Functions

    //- Variables summed over partitions
        //- Total density
        tmp<volScalarField> rho() const;
        
        //- Equivalent potential temperature
        tmp<volScalarField> thetae() const;
        
        //- Density weighted velocity on faces
        tmp<surfaceVectorField> rhoUf() const;
        
        //- Density weighted theta
        tmp<volScalarField> rhoTheta() const;

        //- Update theta and T from Exner in each partition
        void updateThetaT(const volScalarField& Exner);

        // Output
        void write();

        // Input (re-read all fields)
        void readUpdate(const volScalarField& Exner);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
