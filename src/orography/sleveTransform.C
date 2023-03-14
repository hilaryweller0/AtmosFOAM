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

#include "addToRunTimeSelectionTable.H"
#include "dualGradeMountain.H"
#include "sleveTransform.H"

defineTypeNameAndDebug(sleveTransform, 0);
addToRunTimeSelectionTable(terrainFollowingTransform, sleveTransform, dict);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sleveTransform::sleveTransform(const dictionary& dict)
:
H("domainHeight", dimLength, readScalar(dict.lookup("domainHeight"))),
coarseScale(readScalar(dict.lookup("coarseScale"))),
fineScale(readScalar(dict.lookup("fineScale"))),
exponent(dict.lookupOrDefault("exponent", scalar(0))),
m(dualGradeMountain::New(dict.subDict("mountain")))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

point sleveTransform::physicalToComputational(const point& p) const
{
    FatalErrorIn("sleveTransform::physicalToComputational(const point& p)")
        << " not implemented" << exit(FatalError);
    return p;
}

point sleveTransform::computationalToPhysical(const point& p) const
{
	FixedList<scalar,2> s, h, b;
	s[0] = coarseScale;
	s[1] = fineScale;
	const scalar n = exponent;

	scalar z = p.z();

	if (z < H.value())
	{
		h[0] = m->coarseHeightAt(p).value();
		h[1] = m->fineHeightAt(p).value();
		for(int i = 0; i < 2; i++)
		{
			b[i] = Foam::sinh
				   (
					   Foam::pow(H.value()/s[i], n)
					 - Foam::pow(p.z()/s[i], n)
				   )
				   /Foam::sinh(Foam::pow(H.value()/s[i], n));
			z += h[i]*b[i];
		}
	}
    return point(p.x(), p.y(), z);
}

// ************************************************************************* //
