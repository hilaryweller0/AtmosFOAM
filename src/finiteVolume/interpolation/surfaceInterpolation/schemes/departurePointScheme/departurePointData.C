/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "departurePointData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extendedCentredCellToCellStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class ApproxType>
Foam::departurePointData<ApproxType>::departurePointData
(
    const fvMesh& mesh,
    const extendedCentredCellToCellStencil& stencil,
    const surfaceVectorField& depPoints
)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, departurePointData<ApproxType> >(mesh),
    stencil_(stencil),
    departurePoints_(depPoints),
    pointInCell_(mesh.nFaces()),
    coeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing departurePointData<ApproxType>" << endl;
    }

    setPoints();

    if (debug)
    {
        Info<< "departurePointData<ApproxType>::departurePointData() :"
            << "Finished constructing departurePointData"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ApproxType>
void Foam::departurePointData<ApproxType>::setPoints()
{
    const fvMesh& mesh = this->mesh();

//    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
//    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    ApproxType approxType;

    // Diagnostics of the stencil to write out
    surfaceScalarField stencilSizeField
    (
        IOobject("stencilSize", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("s", dimless, scalar(1))
    );
    surfaceScalarField orderField
    (
        IOobject("order", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("s", dimless, scalar(1))
    );
    surfaceScalarField maxStencilWeight
    (
        IOobject("maxStencilWeight", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("s", dimless, scalar(0))
    );
    surfaceScalarField upwindStencilWeight
    (
        IOobject("upwindStencilWeight", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("s", dimless, scalar(0))
    );
    
    // Find the cell containing every departure point
    forAll(departurePoints_, ip)
    {
        const point& p = departurePoints_[ip];
        label celli = mesh.findCell(p);
        if (celli == -1)
        {
            celli = mesh.findNearestCell(p);
        }
        pointInCell_[ip] = celli;
        stencilSizeField[ip] = stencil_.stencil()[celli].size();
    }
    
    // For each departure point, find the nearest face of the cell and
    // the cell on the other side and find the location of that cell
    // in the stencil
    labelList nextCell(departurePoints_.size(), 0);
    forAll(nextCell, ip)
    {
        const label celli = pointInCell_[ip];
        const cell& faces = mesh.cells()[celli];
        scalar minDist = GREAT;
        label minFace = -1;
        forAll(faces, fi)
        {
            const label faceI = faces[fi];
            if (faceI < mesh.nInternalFaces())
            {
                scalar dist = abs
                (
                    (mesh.Cf()[faceI] - departurePoints_[ip])
                  & mesh.Sf()[faceI]
                );
                if (dist < minDist)
                {
                    minDist = dist;
                    minFace = faceI;
                }
            }
        }
        // Now that we have found the closest face, find the 
        // neighbouring cell and find its position in the stencil
        if (minFace >= 0)
        {
            label neiCell = mesh.owner()[minFace] == celli
                          ? mesh.neighbour()[minFace]
                          : mesh.owner()[minFace];

            forAll(stencil_.stencil()[celli], i)
            {
                if (stencil_.stencil()[celli][i] == neiCell)
                {
                    nextCell[ip] = i;
                }
            }
        }
    }

    // Collect the cell centres for all of the stencils
    List<List<point> > Cfld;
    stencil_.collectData(mesh.C(), Cfld);
    
    // Calculate the weights for each stencil
    forAll(departurePoints_, ip)
    {
        label celli = pointInCell_[ip];
        pointField Cs(Cfld[celli]);
        orderField[ip] = approxType.calcWeights
        (
            coeffs_[ip], Cs, departurePoints_[ip], mesh.nSolutionD(),
            nextCell[ip]
        );
        maxStencilWeight[ip] = max(coeffs_[ip]);
        upwindStencilWeight[ip] = coeffs_[ip][0];
    }
    
    // And the boundary values?
    
    // Write out stencil diagnostics
    stencilSizeField.write();
    orderField.write();
    maxStencilWeight.write();
    upwindStencilWeight.write();
}


template<class ApproxType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> > Foam::departurePointData<ApproxType>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    // Collect internal and boundary values
    List<List<Type> > stencilFld;
    stencil_.collectData(vf, stencilFld);
    
    // Declare temporary to return
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(vf.name(), mesh.time().timeName(), mesh),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();
    
    // Sum the contributions for each face from the cell stencils
    forAll(sf, ip)
    {
        label celli = pointInCell_[ip];
        const List<Type>& sfFld = stencilFld[celli];
        const List<scalar>& stWeight = coeffs_[ip];
        
        forAll(sfFld, i)
        {
            sf[ip] += stWeight[i]*sfFld[i];
        }
    }
    
    // Boundary values?
    
    return tsf;
}


template<class ApproxType>
bool Foam::departurePointData<ApproxType>::movePoints()
{
    setPoints();
    return true;
}

// ************************************************************************* //
