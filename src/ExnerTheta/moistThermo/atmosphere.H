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
    A list of fluids with phase changes between them and functions involved in
    combinging the fluids and the equation of state for the mixture

SourceFiles
    atmosphere.C

\*---------------------------------------------------------------------------*/

#ifndef atmosphere_H
#define atmosphere_H

#include "baseAtmosphere.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class atmosphere Declaration
\*---------------------------------------------------------------------------*/

class atmosphere
:
    public baseAtmosphere
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
        //- Volume fraction occupied by dry air
        tmp<volScalarField> volAir() const;
        
        // Pressure from equation of state
        tmp<volScalarField> pFromT(const volScalarField& T) const;
        
        // Exner from equation of state
        tmp<volScalarField> ExnerFromTheta(const volScalarField& theta) const;
        
        // Exner from equation of state
        void setExnerFromTheta
        (
            volScalarField& Exner,
            const volScalarField& theta
        ) const;
        
        // rho from eqn of state
        tmp<volScalarField> rhoFromP
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;
        
        // rho from eqn of state
        tmp<volScalarField> rhoFromExner
        (
            const volScalarField& Exner,
            const volScalarField& theta
        ) const;
        
        // Source term in the theta equation (divided by rhoa theta)
        tmp<volScalarField> thetaSource
        (
            const volScalarField& T,
            const volScalarField& divu
        ) const;
        
        tmp<volScalarField> thetae(const volScalarField& T) const;
        
        //- Calculate temperautre, T, from uniform equivalent potential
        //- temperature and uniform total water mixing ratio, rt
        void TfromThetae
        (
            volScalarField& T,
            const dimensionedScalar& thetae0,
            const scalar rt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
