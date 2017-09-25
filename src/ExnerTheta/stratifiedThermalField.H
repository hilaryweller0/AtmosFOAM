/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::stratifiedThermalField

Description
    A stratified potential temperature field defined as a tracerField
    implementation.  The stratification is defined by the Brunt–Väisälä
    frequency configured in environmentalProperties and the surface
    temperature configured in thermophysicalProperties.

SourceFiles
    stratifiedThermalField.C

\*---------------------------------------------------------------------------*/

#ifndef stratifiedThermalField_H
#define stratifiedThermalField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"

namespace Foam
{

class stratifiedThermalField;

/*---------------------------------------------------------------------------*\
                         Class stratifiedThermalField Declaration
\*---------------------------------------------------------------------------*/

class stratifiedThermalField
{
    // Private data

        const dimensionedVector g;
        const dimensionedScalar T0;
        const scalarList nLayers;
        const scalarList zLayers;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        stratifiedThermalField(const stratifiedThermalField&);

        //- Disallow default bitwise assignment
        void operator=(const stratifiedThermalField&);


public:

    // Constructors

        stratifiedThermalField(
                const IOdictionary& environmentalProperties,
                const dimensionedVector g,
                const dimensionedScalar T0);

    // Selectors

    // Member Functions

        scalar tracerAt(const point& p, const Time& t) const;
        vector gradAt(const point& p, const Time& t) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
