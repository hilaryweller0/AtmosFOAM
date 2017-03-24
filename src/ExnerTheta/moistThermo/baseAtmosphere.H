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
    Foam::baseAtmosphere

Description
    A list of fluids with phase changes between them and functions involved in
    combinging the fluids and the equation of state for the mixture

SourceFiles
    baseAtmosphere.C

\*---------------------------------------------------------------------------*/

#ifndef baseAtmosphere_H
#define baseAtmosphere_H

#include "fluidSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class baseAtmosphere Declaration
\*---------------------------------------------------------------------------*/

class baseAtmosphere
:
    public PtrList<fluidSpecie>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        baseAtmosphere(const baseAtmosphere&);

        //- Disallow default bitwise assignment
        void operator=(const baseAtmosphere&);


public:

    // Constructors

        //- Construct
        baseAtmosphere
        (
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict
        );


    //- Destructor
    ~baseAtmosphere();


    // Member Functions

        // Access
        
        //- Volume fraction occupied by gases
        tmp<volScalarField> volGas() const;
        
        //- Total density
        tmp<volScalarField> sumDensity() const;
        
        //- Sum of partial pressures
        tmp<volScalarField> sumPressure(const volScalarField& T) const;
        
        //- Density weighted total gas constant
        tmp<volScalarField> rhoR() const;
        
        //- Density weighted total heat capacity (constant pressure)
        tmp<volScalarField> rhoCp() const;
        
        //- Density weighted total heat capacity (constant volume)
        tmp<volScalarField> rhoCv() const;
        
        // Output
        void write();

        // Input (re-read densities and liquid fractions)
        void readUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
