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

Class
    Foam::linearCP

Description
    Linear interpolation using a least squares fit of a_1 + a_2 z to
    surrounding faces.  Vertical face values are not included in the fit.

SourceFiles
    linearCP.C

\*---------------------------------------------------------------------------*/

#ifndef linearCP_H
#define linearCP_H

#include "volInterpolationScheme.H"
#include "volFields.H"
#include "singleFaceToCellStencil.H"
#include "extendedCentredFaceToCellStencil.H"
#include "gravity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearCP Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearCP
:
    public volInterpolationScheme<Type>
{
public:
    TypeName("linearCP");

    linearCP(const fvMesh& mesh)
    :
        volInterpolationScheme<Type>(mesh),
        f2cStencil(mesh),
        stencilDescription(f2cStencil),
        coeffs(mesh.nCells()),
        g(mesh)
    {
        initCoeffs(mesh);
    }

    linearCP(const fvMesh& mesh, Istream&)
    :
        volInterpolationScheme<Type>(mesh),
        f2cStencil(mesh),
        stencilDescription(f2cStencil),
        coeffs(mesh.nCells()),
        g(mesh)
    {
        initCoeffs(mesh);
    }

    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>
    interpolate
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const;

private:
    linearCP(const linearCP&);
    void operator=(const linearCP&);

    void initCoeffs(const fvMesh& mesh);
    void calculateInterpolationCoeffs
    (
        const point& origin,
        const labelList& stencilFaceIndices,
        const List<point>& stencil,
        scalarList& coeffs
    );

    const singleFaceToCellStencil f2cStencil;
    const extendedCentredFaceToCellStencil stencilDescription;
    List<scalarList> coeffs;
    const gravity g;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "linearCP.C"
#endif

#endif

// ************************************************************************* //
