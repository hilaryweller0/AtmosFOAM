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
    Foam::dualGradeMountain

Description

SourceFiles
    dualGradeMountain.C

\*---------------------------------------------------------------------------*/

#ifndef dualGradeMountain_H
#define dualGradeMountain_H

#include "mountain.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class dualGradeMountain Declaration
\*---------------------------------------------------------------------------*/

class dualGradeMountain
:
    public mountain
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //dualGradeMountain(const dualGradeMountain&);

        //- Disallow default bitwise assignment
        //void operator=(const dualGradeMountain&);


public:

    // Selectors

        declareRunTimeSelectionTable
        (
            autoPtr,
            dualGradeMountain,
            dict,
            (
                const dictionary& dict
            ),
            (dict)
        );

        static autoPtr<dualGradeMountain> New(const dictionary& dict);

    // Member Functions

        dimensionedScalar heightAt(const point& p) const;

        virtual dimensionedScalar coarseHeightAt(const point& p) const = 0;

        virtual dimensionedScalar fineHeightAt(const point& p) const = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
