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

#include "fvc.H"
#include "fcfBilinearFit.H"
#include "SVD.H"

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvsPatchField,
        Foam::surfaceMesh
    >
> Foam::fv::fcfBilinearFit<Type>::operator()
(
    const Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>& s
) const
{
    return stencilDescription.weightedSum(s, coeffs);
}

template<class Type>
void Foam::fv::fcfBilinearFit<Type>::initCoeffs
(
    const fvMesh& mesh
)
{
    List<List<point> > stencilPoints(mesh.nFaces());

    stencilDescription.collectData
    (
        mesh.Cf(),
        stencilPoints
    );

    forAll(stencilPoints, stencilForFaceI)
    {
        const List<point>& stencil = stencilPoints[stencilForFaceI];

        if (stencil.size() >= 3)
        {
            calculateGradCoeffs
            (
                stencilDescription.elements()[stencilForFaceI],
                stencil,
                coeffs[stencilForFaceI]
            );
        }
        else
        {
            // not enough faces to fit phi = a_1 + a_2 x + a_3 z
            Info << "ignoring faceI " << stencilForFaceI << " stencil of size " << stencil.size() << endl;
            forAll(stencil, faceI)
            {
                coeffs[stencilForFaceI].append(vector(0, 0, 0));
            }
        }
    }
}

template<class Type>
void Foam::fv::fcfBilinearFit<Type>::calculateGradCoeffs
(
    const labelList& stencilFaceIndices,
    const List<point>& stencil,
    List<vector>& coeffs
)
{
    label nonVerticalFaces = 0;
    boolList include(stencilFaceIndices.size());
    scalarList multipliers(stencil.size());
    forAll(stencilFaceIndices, i)
    {
        multipliers[i] = mag(g.unitFaceNormal()[stencilFaceIndices[i]]);
        include[i] = mag(g.unitFaceNormal()[stencilFaceIndices[i]]) > 1e-12;
        if (include[i]) nonVerticalFaces++;
    }

    scalarRectangularMatrix B(nonVerticalFaces, 3);

    label n = 0;
    const point& origin = stencil.first();
    forAll(stencil, i)
    {
        if (include[i])
        {
            const point& p = stencil[i] - origin;

            B(n, 0) = multipliers[i] * 1;
            B(n, 1) = multipliers[i] * p.x();
            B(n, 2) = multipliers[i] * p.z();

            n++;
        }
    }

    const SVD svd(B);
    const scalarRectangularMatrix& Binv = svd.VSinvUt();

    n = 0;
    forAll(stencil, i)
    {
        if (include[i])
        {            
            coeffs.append(vector(multipliers[i]*Binv(1, n), 0, multipliers[i]*Binv(2, n)));
            n++;
        }
        else
        {
            coeffs.append(vector(0, 0, 0));
        }
    }
}
