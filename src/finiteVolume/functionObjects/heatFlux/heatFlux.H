/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::functionObjects::heatFlux

Description
    Calculates and write the heat-flux as the volScalarField field 'heatFlux'.

SourceFiles
    heatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_heatFlux_H
#define functionObjects_heatFlux_H

#include "fieldsExpression.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class heatFlux Declaration
\*---------------------------------------------------------------------------*/

class heatFlux
:
    public fieldsExpression
{
    // Private Data

    //- Names of the velocity, flux, density and temperature fields
    word UName_, fluxName_, rhoName_, TName_;


    // Private Member Functions

        //- Calculate the heat flux
        virtual bool calc();
public:

    //- Runtime type information
    TypeName("heatFlux");


    // Constructors

        //- Construct from Time and dictionary
        heatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~heatFlux();


    // Member Functions

        //- Read the heatFlux data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
