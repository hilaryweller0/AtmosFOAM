/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "fvcDet.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> det(const volTensorField& A)
{
    const fvMesh& mesh = A.mesh();

    tmp<volScalarField> tdet
    (
        new volScalarField
        (
            IOobject("det("+A.name()+")", A.instance(), mesh),
            mesh,
            dimensionedScalar("", pow3(A.dimensions()), scalar(0))
        )
    );
    
    volScalarField& d = tdet.ref();
    
    forAll(d, cellI)
    {
        d[cellI] = det(A[cellI]);
    }
    forAll(mesh.boundary(), patchi)
    {
        forAll(d.boundaryField()[patchi], faci)
        {
            d.boundaryFieldRef()[patchi][faci]
                 = det(A.boundaryField()[patchi][faci]);
        }
    }
    
    return tdet;
}


tmp<volScalarField> det(const tmp<volTensorField>& tA)
{
    tmp<volScalarField> tdet
    (
        fvc::det(tA())
    );
    tA.clear();
    return tdet;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
