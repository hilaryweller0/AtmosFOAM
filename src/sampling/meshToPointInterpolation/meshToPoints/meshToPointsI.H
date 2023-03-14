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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class approxType>
template<class Type>
inline Field<Type> Foam::meshToPoints<approxType>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    Field<Type> rVal(mesh_.nPoints(), pTraits<Type>::zero);
    forAll(rVal, ip)
    {
        const labelList& pointCells = mesh_.pointCells()[ip];
        forAll(pointCells, ic)
        {
            label cellI = pointCells[ic];
            rVal[ip] += weights_[ip][ic]*vf[cellI];
        }
    }
    
    return rVal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
