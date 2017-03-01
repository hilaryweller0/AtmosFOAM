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

Class
    Foam::fcfBilinearFit

Description
    A surface grad scheme that uses a face-cell-face stencil and a bilinear
    least squares fit.

SourceFiles
    fcfBilinearFit.C

\*---------------------------------------------------------------------------*/

#ifndef fcfBilinearFit_H
#define fcfBilinearFit_H

#include "sGradScheme.H"
#include "extendedCentredFaceToFaceStencil.H"
#include "CFCFaceToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

template<class Type>
class fcfBilinearFit
:
    public fv::sGradScheme<Type>
{
public:
    TypeName("fcfBilinearFit");

    fcfBilinearFit(const fvMesh& mesh)
    :
        sGradScheme<Type>(mesh),
        coeffs(mesh.nFaces()),
        f2fStencil(mesh),
        stencil(f2fStencil)
    {
        initCoeffs(mesh);
    }

    fcfBilinearFit(const fvMesh& mesh, Istream&)
    :
        sGradScheme<Type>(mesh),
        coeffs(mesh.nFaces()),
        f2fStencil(mesh),
        stencil(f2fStencil)
    {
        initCoeffs(mesh);
    }

    virtual tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvsPatchField, surfaceMesh>
    > operator()
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const;
    
private:
    fcfBilinearFit(const fcfBilinearFit&);
    void operator=(const fcfBilinearFit&);

    void initCoeffs(const fvMesh& mesh);

    List<List<vector > > coeffs;
    const CFCFaceToFaceStencil f2fStencil;
    const extendedCentredFaceToFaceStencil stencil;
};

}

}

#ifdef NoRepository
    #include "fcfBilinearFit.C"
#endif

#endif
