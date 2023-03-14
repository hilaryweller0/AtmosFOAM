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

#include "meshToPointField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class ApproxType, class ExtendedStencil>
Foam::meshToPointField<Type, ApproxType, ExtendedStencil>::meshToPointField
(
    const fvMesh& mesh,
    const pointField& pts,
    const ExtendedStencil& cellStencil__
)
:
    mesh_(mesh),
    cellStencil_(cellStencil__),
    pointStencil_(),
    weights_()
{
    setPoints(pts);
}


template<class Type, class ApproxType, class ExtendedStencil>
Foam::meshToPointField<Type, ApproxType, ExtendedStencil>::meshToPointField
(
    const fvMesh& mesh,
    const ExtendedStencil& cellStencil__
)
:
    mesh_(mesh),
    cellStencil_(cellStencil__),
    pointStencil_(),
    weights_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class ApproxType, class ExtendedStencil>
Foam::meshToPointField<Type, ApproxType, ExtendedStencil>::~meshToPointField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type, class ApproxType, class ExtendedStencil>
void Foam::meshToPointField
<
    Type, ApproxType, ExtendedStencil
>::setPoints(const pointField& pts)
{
    weights_.setSize(pts.size());
    pointStencil_.setSize(pts.size());
    
    ApproxType approxType;

    // Write out the stencil size and order of accuracy for every face
    surfaceScalarField stencilSizeField
    (
        IOobject("stencilSize", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );
    surfaceScalarField orderField
    (
        IOobject("order", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );

    // Set the stencil for each point depending on the cell it is in
    forAll(pts, ip)
    {
        const point& p = pts[ip];
        label celli = mesh_.findCell(p);
        if (celli == -1)
        {
            celli = mesh_.findNearestCell(p);
        }
        pointStencil_[ip] = cellStencil_.stencil()[celli];
        stencilSizeField[ip] = pointStencil_[ip].size();
    }

    // Collect internal and boundary values
    List<List<point> > Cfld;
    extendedCellToFaceStencil::collectData
    (
        cellStencil().map(), pointStencil_, mesh_.C(), Cfld
    );

    // Calculate the weights for each stencil
    forAll(pts, ip)
    {
        pointField Cs(Cfld[ip]);
        orderField[ip] = approxType.calcWeights
        (
            weights_[ip], Cs, pts[ip], mesh_.nSolutionD()
        );
        
//        if (ip == 438)
//        {
//            Info << "pt = " << pts[ip] << " Cs = " << Cs << endl;
//        }
    }
//    label ip = 438;
//    Info << "Face " << ip << " weights " << weights_[ip] << endl;
    stencilSizeField.write();
    orderField.write();
}

// ************************************************************************* //
