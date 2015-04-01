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

#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class approxType>
template<class Type>
inline Field<Type> Foam::faceToPoint<approxType>::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    Field<Type> rVal(mesh_.nPoints(), pTraits<Type>::zero);
    labelList pointFaces(3);
    forAll(rVal, ip)
    {
        //const labelList& pointFaces = mesh_.pointFaces()[ip];
        pointFaces[0] = mesh_.pointFaces()[ip][0];
        pointFaces[1] = mesh_.pointFaces()[ip][1];
        pointFaces[2] = mesh_.pointFaces()[ip][2];
        forAll(pointFaces, iff)
        {
            label faceI = pointFaces[iff];
            if (faceI < mesh_.nInternalFaces())
            {
                rVal[ip] += weights_[ip][iff]*sf[faceI];
            }
        }
    }
    
    return rVal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
