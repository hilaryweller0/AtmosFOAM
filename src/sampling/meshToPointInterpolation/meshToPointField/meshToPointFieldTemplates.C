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

#include "extendedCellToCellStencil.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class ApproxType, class ExtendedStencil>
Foam::tmp<Foam::Field<Type> >
Foam::meshToPointField<Type, ApproxType, ExtendedStencil>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (weights_.size() == 0)
    {
        FatalErrorIn
        (
            "meshToPointField<ApproxType, ExtendedStencil>::interpolate"
        )
            << " points to interpolate onto not set"
            << exit(FatalError);
    }

    // Collect internal and boundary values
    List<List<Type> > stencilFld;
    extendedCellToFaceStencil::collectData
    (
        cellStencil().map(), pointStencil_, vf, stencilFld
    );

    // Declare temporary field to return
    tmp<Field<Type> > tpf(new Field<Type>(pointStencil_.size(), pTraits<Type>::zero));
    Field<Type> & pf = tpf();
    
    forAll(pf, ip)
    {
        const List<Type>& stField = stencilFld[ip];
        const scalarList& stWeight = weights_[ip];

        forAll(stField, i)
        {
            pf[ip] += stWeight[i]*stField[i];
        }
    }
    
    return tpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
