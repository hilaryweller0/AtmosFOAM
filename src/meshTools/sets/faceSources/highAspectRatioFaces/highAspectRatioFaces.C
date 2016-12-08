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

#include "highAspectRatioFaces.H"
#include "polyMesh.H"
#include "faceSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(highAspectRatioFaces, 0);

addToRunTimeSelectionTable(topoSetSource, highAspectRatioFaces, word);

addToRunTimeSelectionTable(topoSetSource, highAspectRatioFaces, istream);

}


Foam::topoSetSource::addToUsageTable Foam::highAspectRatioFaces::usage_
(
    highAspectRatioFaces::typeName,
    "\n    Usage: highAspectRatioFaces\n\n"
    "    Select faces whose removal would lead to a cell with lower\n"
    "    aspect ratio than the original cells either side of the face\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::highAspectRatioFaces::highAspectRatioFaces
(
    const polyMesh& mesh
)
:
    topoSetSource(mesh)
{

}


// Construct from dictionary
Foam::highAspectRatioFaces::highAspectRatioFaces(const polyMesh& mesh, const dictionary& dict)
:
    topoSetSource(mesh)
{}


// Construct from Istream
Foam::highAspectRatioFaces::highAspectRatioFaces(const polyMesh& mesh, Istream& is)
:
    topoSetSource(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highAspectRatioFaces::~highAspectRatioFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::highAspectRatioFaces::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding high aspect ratio faces "
            << endl;

        forAll(mesh_.faceNeighbour(), facei)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            scalar Vown = mesh_.cellVolumes()[own];
            scalar Vnei = mesh_.cellVolumes()[nei];
            scalar magSf = mag(mesh_.faceAreas()[facei]);
            
            if (pow(Vown + Vnei, 2/3.) < magSf)
            {
                set.insert(facei);
            }
        }
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing high aspect ratio faces "
            << endl;


        DynamicList<label> toBeRemoved(set.size()/10);

        forAllConstIter(topoSet, set, iter)
        {
            const label facei = iter.key();

            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            scalar Vown = mesh_.cellVolumes()[own];
            scalar Vnei = mesh_.cellVolumes()[nei];
            scalar magSf = mag(mesh_.faceAreas()[facei]);

            if (pow(Vown + Vnei, 2/3.) < magSf)
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
