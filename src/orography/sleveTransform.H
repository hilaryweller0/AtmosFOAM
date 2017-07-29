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
    Foam::sleveTransform

Description
    Applies the Smooth LEvel VErtical (SLEVE) transform as specified by
    https://doi.org/dtdbcm, Schaer et al 2002

    Requires a dualGradeMountain implementation in order to separate
    coarse and fine-scale topography.

SourceFiles
    sleveTransform.C

\*---------------------------------------------------------------------------*/

#ifndef sleveTransform_H
#define sleveTransform_H

#include "dualGradeMountain.H"
#include "terrainFollowingTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sleveTransform Declaration
\*---------------------------------------------------------------------------*/

class sleveTransform final
:
    public terrainFollowingTransform
{
    // Private data

    const dimensionedScalar H;
    const scalar coarseScale;
    const scalar fineScale;
    const scalar exponent;
    const autoPtr<dualGradeMountain> m;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sleveTransform(const sleveTransform&);

        //- Disallow default bitwise assignment
        void operator=(const sleveTransform&);


public:

    TypeName("sleve");

    // Constructors

        sleveTransform(const dictionary& dict);

    // Member Functions

        point physicalToComputational(const point& p) const;

        point computationalToPhysical(const point& p) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
