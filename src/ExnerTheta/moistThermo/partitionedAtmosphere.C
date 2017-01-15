/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "partitionedAtmosphere.H"
#include "moreListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::partitionedAtmosphere
(
    const wordList& partitionNames,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    PtrList<partition>(partitionNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new partition
            (
                partitionNames[ip],
                partNames,
                mesh,
                dict
            )
        );
    }
}


Foam::partitionedAtmosphere::partitionedAtmosphere
(
    const wordList& partitionNames,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict,
    const volScalarField& Exner
)
:
    PtrList<partition>(partitionNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new partition
            (
                partitionNames[ip],
                partNames,
                mesh,
                dict,
                Exner
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partitionedAtmosphere::~partitionedAtmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::partitionedAtmosphere::updateThetaT(const volScalarField& Exner)
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.updateThetaT(Exner);
    }
}


void Foam::partitionedAtmosphere::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.write();
    }
}


void Foam::partitionedAtmosphere::readUpdate(const volScalarField& Exner)
{
    for(label ip = 0; ip < size(); ip++)
    {
        partition& part = operator[](ip);
        part.readUpdate(Exner);
    }
}


// ************************************************************************* //
