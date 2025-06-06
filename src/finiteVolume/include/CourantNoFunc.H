/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    CourantNoFunc

Description
    Calculates and returns a volScalarField of Courant Numbers.

\*---------------------------------------------------------------------------*/

#ifndef CourantNoFunc_H
#define CourantNoFunc_H

namespace Foam
{
inline tmp<volScalarField> CourantNo
(
    const surfaceScalarField& phi,
    const dimensionedScalar& dt
)
{
    const fvMesh& mesh = phi.mesh();

    tmp<volScalarField> tCourantNo
    (
        new volScalarField
        (
            IOobject("Co", mesh.time().name(), mesh),
            mesh,
            dimensionedScalar(dimless, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
    tCourantNo->primitiveFieldRef() = (sumPhi/mesh.V().primitiveField())
                                    *(0.5*dt.value());
    tCourantNo.ref().correctBoundaryConditions();

    return tCourantNo;
}
}
#endif

