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

#include "fvcCof.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> cof(const volTensorField& A)
{
    const fvMesh& mesh = A.mesh();

    tmp<volTensorField> tcof
    (
        new volTensorField
        (
            IOobject("cof("+A.name()+")", A.instance(), mesh),
            A
        )
    );
    
    volTensorField& C = tcof.ref();
    
    forAll(A, cellI)
    {
        C[cellI] = cof(A[cellI]);
    }
    forAll(mesh.boundary(), patchi)
    {
        forAll(A.boundaryField()[patchi], faci)
        {
            C.boundaryFieldRef()[patchi][faci]
                 = cof(A.boundaryField()[patchi][faci]);
        }
    }
    
    return tcof;
}


tmp<volTensorField> cof(const tmp<volTensorField>& tA)
{
    tmp<volTensorField> tcof
    (
        fvc::cof(tA())
    );
    tA.clear();
    return tcof;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
