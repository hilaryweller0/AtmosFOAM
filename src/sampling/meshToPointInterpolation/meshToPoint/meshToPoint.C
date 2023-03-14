/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "meshToPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPoint<approxType>::meshToPoint
(
    const point& p, const fvMesh& mesh
)
:
    point(p),
    mesh_(mesh),
    stencil_(),
    weights_()
{
    label celli = mesh_.findNearestCell(p);
    label stencilSize = mesh_.cellCells()[celli].size()+1;
    stencil_.setSize(stencilSize);
    weights_.setSize(stencilSize);
    
    // Set the stencil to be the central cell plus its neighbours
    stencil_[0] = celli;
    for(label is = 1; is < stencilSize; is++)
    {
        stencil_[is] = mesh_.cellCells()[celli][is-1];
    }
    
    // Create the weights (depending on the approxType)
    approxType thisApproxType;
    thisApproxType.calcWeights(weights_, mesh_, stencil_, p);
}


template<class approxType>
Foam::meshToPoint<approxType>::meshToPoint(const fvMesh& mesh)
:
    point(vector::zero),
    mesh_(mesh),
    stencil_(),
    weights_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPoint<approxType>::~meshToPoint()
{}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class approxType>
void Foam::meshToPoint<approxType>::setPoint(const point& pt)
{
    this->point = pt;
    
    label celli = mesh_.findNearestCell(pt);
    label stencilSize = mesh_.cellCells()[celli].size()+1;
    stencil_.setSize(stencilSize);
    weights_.setSize(stencilSize);
    
    // Set the stencil to be the central cell plus its neighbours
    stencil_[0] = celli;
    for(label is = 1; is < stencilSize; is++)
    {
        stencil_[is] = mesh_.cellCells()[celli][is-1];
    }
    
    // Create the weights (depending on the approxType)
    approxType thisApproxType;
    thisApproxType.calcWeights(weights_, mesh_, stencil_, pt);
}

// ************************************************************************* //
