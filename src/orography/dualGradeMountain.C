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

\*---------------------------------------------------------------------------*/

#include "dualGradeMountain.H"

namespace Foam
{
defineRunTimeSelectionTable(dualGradeMountain, dict);

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dualGradeMountain> dualGradeMountain::New(const dictionary& dict)
{
    const word mountainType(dict.lookup("type"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(mountainType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dualGradeMountain::New(const dictionary&)"
        ) << "Unknown type "
          << mountainType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<dualGradeMountain>
    (
       cstrIter()(dict)
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

dimensionedScalar dualGradeMountain::heightAt(const point& p) const
{
    return coarseHeightAt(p) + fineHeightAt(p);
}

// ************************************************************************* //
}
