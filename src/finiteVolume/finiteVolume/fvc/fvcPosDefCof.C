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

#include "fvcPosDefCof.H"
#include "fvcCof.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> posDefCof
(
    const volTensorField& A,
    const scalar minEig
)
{
    const fvMesh& mesh = A.mesh();

    tmp<volTensorField> tcof
    (
        new volTensorField
        (
            IOobject("cof("+A.name()+")", A.instance(), mesh),
            fvc::cof(A)
        )
    );
    
    volTensorField& C = tcof.ref();
    
    // Ensure that all eigenvalues are positive
    
    forAll(C, cellI)
    {
        scalar smallestEig = eigenValues(C[cellI])[0];
        if (smallestEig < minEig)
        {
            C[cellI] += (minEig - smallestEig)*tensor::I;
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        forAll(C.boundaryField()[patchi], faci)
        {
            scalar smallestEig
                 = eigenValues(C.boundaryFieldRef()[patchi][faci])[0];
            C.boundaryFieldRef()[patchi][faci]
                 += (minEig - smallestEig)*tensor::I;
        }
    }
    
    return tcof;
}


tmp<volTensorField> posDefCof
(
    const tmp<volTensorField>& tA,
    const scalar minEig
)
{
    tmp<volTensorField> tcof
    (
        fvc::posDefCof(tA(), minEig)
    );
    tA.clear();
    return tcof;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
