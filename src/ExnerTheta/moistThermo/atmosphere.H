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
    Foam::atmosphere

Description
    A list of fluids with phase changes between them and fucntions involved in
    combinging the fluids and the equation of state for the mixture

SourceFiles
    atmosphere.C

\*---------------------------------------------------------------------------*/

#ifndef atmosphere_H
#define atmosphere_H

#include "fluidSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class atmosphere Declaration
\*---------------------------------------------------------------------------*/

class atmosphere
:
    public PtrList<fluidSpecie>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        atmosphere(const atmosphere&);

        //- Disallow default bitwise assignment
        void operator=(const atmosphere&);


public:

    // Constructors

        //- Construct
        atmosphere
        (
            const wordList& partNames,
            const fvMesh& mesh,
            const dictionary dict
        );


    //- Destructor
    ~atmosphere();


    // Member Functions

        // Access
        //- Volume fraction occupied by gases
        tmp<volScalarField> volGas();
        //- Volume fraction occupied by dry air
        tmp<volScalarField> volAir();
        //- Total density
        tmp<volScalarField> sumDensity();
        //- Density weighted total gas constant
        tmp<volScalarField> rhoR();
        //- Density weighted total heat capacity (constant pressure)
        tmp<volScalarField> rhoCp();
        //- Density weighted total heat capacity (constant volume)
        tmp<volScalarField> rhoCv();
        //- Sum of partial pressures
        tmp<volScalarField> sumPressure(const volScalarField& T);
        // Pressure from equation of state
        tmp<volScalarField> pFromT(const volScalarField& T);
        // Exner from equation of state
        tmp<volScalarField> ExnerFromTheta(const volScalarField& theta);
        // rho from eqn of state
        tmp<volScalarField> rhoFromP
        (
            const volScalarField& p,
            const volScalarField& T
        );
        // Source term in the theta equation (divided by rhoa theta)
        tmp<volScalarField> thetaSource
        (
            const volScalarField& T,
            const volScalarField& divu
        );
        
        tmp<volScalarField> thetae(const volScalarField& T);
        
        //- Calculate temperautre, T, from uniform equivalent potential
        //- temperature and uniform total water mixing ratio, rt
        void TfromThetae
        (
            volScalarField& T,
            const dimensionedScalar& thetae0,
            const scalar rt
        );

        // Output
        void write();

        // Input (re-read densities and liquid fractions
        void readUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
