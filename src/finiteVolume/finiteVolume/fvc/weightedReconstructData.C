/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "weightedReconstructData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(weightedReconstructData, 0);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::weightedReconstructData::weightedReconstructData
(
    const fvMesh& mesh,
    const scalar boundaryWeight__
)
:
    DemandDrivenMeshObject<fvMesh, Foam::MoveableMeshObject, weightedReconstructData>(mesh),
    boundaryWeight_(boundaryWeight__),
    faceWeights_
    (
        IOobject("faceWeights", this->mesh().pointsInstance(), this->mesh()),
        this->mesh(),
        dimensionedVector(dimless, Zero)
    ),
    invTensor_
    (
        IOobject("invTensor", this->mesh().pointsInstance(), this->mesh()),
        this->mesh(),
        dimensionedTensor(dimless/dimArea, Zero)
    )
{
    calcData();
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::weightedReconstructData::~weightedReconstructData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::weightedReconstructData::calcData()
{
    if (debug)
    {
        InfoInFunction << "Calculating reconstruction daa" << endl;
    }

    faceWeights_ = this->mesh().Sf()/this->mesh().magSf();
    forAll(faceWeights_.boundaryField(), patchi)
    {
        if (!this->mesh().boundaryMesh()[patchi].coupled())
        {
            faceWeights_.boundaryFieldRef()[patchi] *= boundaryWeight_;
        }
    }
    
    invTensor_ = inv(fvc::surfaceSum(faceWeights_*this->mesh().Sf()));

    if (debug)
    {
        InfoInFunction
            << "Finished calculating reconstruction data" << endl;
    }
}


bool Foam::weightedReconstructData::movePoints()
{
    calcData();
    return true;
}

// ************************************************************************* //
