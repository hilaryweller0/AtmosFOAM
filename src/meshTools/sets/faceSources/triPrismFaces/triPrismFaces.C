/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "triPrismFaces.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "prismMatcher.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triPrismFaces, 0);

addToRunTimeSelectionTable(topoSetSource, triPrismFaces, word);

addToRunTimeSelectionTable(topoSetSource, triPrismFaces, istream);

}


Foam::topoSetSource::addToUsageTable Foam::triPrismFaces::usage_
(
    triPrismFaces::typeName,
    "\n    Usage: triPrismFaces\n\n"
    "    Select faces whose owner and neighbour are triangular prisms\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::triPrismFaces::triPrismFaces
(
    const polyMesh& mesh
)
:
    topoSetSource(mesh)
{

}


// Construct from dictionary
Foam::triPrismFaces::triPrismFaces(const polyMesh& mesh, const dictionary& dict)
:
    topoSetSource(mesh)
{}


// Construct from Istream
Foam::triPrismFaces::triPrismFaces(const polyMesh& mesh, Istream& is)
:
    topoSetSource(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triPrismFaces::~triPrismFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triPrismFaces::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    prismMatcher prism;

    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces with owner and neighbour triangular prisms "
            << endl;

        forAll(mesh_.faceNeighbour(), facei)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            const cell& cellOwn = mesh_.cells()[own];
            const cell& cellNei = mesh_.cells()[nei];

            if
            (
                prism.isA(mesh_, own) && prism.isA(mesh_, nei)
             && cellOwn.nFaces() == 5 && cellNei.nFaces() == 5
            )
            {
                set.insert(facei);
            }
        }
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces with owner and neighbour triangular prisms "
            << endl;


        DynamicList<label> toBeRemoved(set.size()/10);

        forAllConstIter(topoSet, set, iter)
        {
            const label facei = iter.key();

            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            const cell& cellOwn = mesh_.cells()[own];
            const cell& cellNei = mesh_.cells()[nei];

            if
            (
                prism.isA(mesh_, own) && prism.isA(mesh_, nei)
             && cellOwn.nFaces() == 5 && cellNei.nFaces() == 5
            )
            {
                toBeRemoved.append(facei);
            }
        }

        forAll(toBeRemoved, i)
        {
            set.erase(toBeRemoved[i]);
        }
    }
}


// ************************************************************************* //
