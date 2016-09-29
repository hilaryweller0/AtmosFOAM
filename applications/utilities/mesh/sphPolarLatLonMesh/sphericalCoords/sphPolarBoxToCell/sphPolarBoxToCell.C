/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "sphPolarBoxToCell.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sphPolarBoxToCell, 0);

addToRunTimeSelectionTable(topoSetSource, sphPolarBoxToCell, word);

addToRunTimeSelectionTable(topoSetSource, sphPolarBoxToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::sphPolarBoxToCell::usage_
(
    sphPolarBoxToCell::typeName,
    "\n    Usage: sphPolarBoxToCell (minLon minLat minR) (maxLon maxLat maxR)\n\n"
    "    Select all cells with cellCentre within lon/lat/radius bounding box\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphPolarBoxToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, cellI)
    {
        polarPoint cp = convertToPolar
        (
            ctrs[cellI], max(scalar(360), TRcorner_.lon())
        );
        if
        (
            cp.lon() >= BLcorner_.lon() && cp.lat() >= BLcorner_.lat() && cp.r() >= BLcorner_.r()
         && cp.lon() <= TRcorner_.lon() && cp.lat() <= TRcorner_.lat() && cp.r() <= TRcorner_.r()
        )
        {
            addOrDelete(set, cellI, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::sphPolarBoxToCell::sphPolarBoxToCell
(
    const polyMesh& mesh,
    const polarPoint& BLcorner,
    const polarPoint& TRcorner
)
:
    topoSetSource(mesh),
    BLcorner_(BLcorner),
    TRcorner_(TRcorner)
{
    if (TRcorner_.lon() < BLcorner_.lon())
    {
        TRcorner_.lon() += 360.;
    }
}


// Construct from dictionary
Foam::sphPolarBoxToCell::sphPolarBoxToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    BLcorner_(dict.lookup("BLcorner")),
    TRcorner_(dict.lookup("TRcorner"))
{
    if (TRcorner_.lon() < BLcorner_.lon())
    {
        TRcorner_.lon() += 360.;
    }
}


// Construct from Istream
Foam::sphPolarBoxToCell::sphPolarBoxToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    BLcorner_(checkIs(is)),
    TRcorner_(checkIs(is))
{
    if (TRcorner_.lon() < BLcorner_.lon())
    {
        TRcorner_.lon() += 360.;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sphPolarBoxToCell::~sphPolarBoxToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphPolarBoxToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with center within spherical polar box " << BLcorner_
            << ' ' << TRcorner_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with center within spherical polar box " << BLcorner_
            << ' ' << TRcorner_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
