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
    Foam::functionObjects::wFlux

Description
    Calculates and write the w-flux as the volScalarField field 'wFlux'.

SourceFiles
    wFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wFlux_H
#define functionObjects_wFlux_H

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
                        Class wFlux Declaration
\*---------------------------------------------------------------------------*/

class wFlux
:
    public fieldsExpression
{
    // Private Data

    //- Names of the velocity, flux, density and temperature fields
    word UName_, fluxName_, rhoName_, pName_;


    // Private Member Functions

        //- Calculate the heat flux
        virtual bool calc();
public:

    //- Runtime type information
    TypeName("wFlux");


    // Constructors

        //- Construct from Time and dictionary
        wFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wFlux();


    // Member Functions

        //- Read the wFlux data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
