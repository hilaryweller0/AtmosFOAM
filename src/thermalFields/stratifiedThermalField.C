/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "stratifiedThermalField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(stratifiedThermalField, 0);
addToRunTimeSelectionTable(tracerField, stratifiedThermalField, dict);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::stratifiedThermalField::checkLayerSizes()
{
    if (nLayers.size()+1 != zLayers.size())
    {
        FatalErrorIn("stratifiedThermalField")
            << " size of BruntVaisalaFrequencies should be"
            << " one smaller than the size of verticalLayers"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stratifiedThermalField::stratifiedThermalField
(
        const dictionary& dict,
        const advectable& velocityField
)
:
    tracerField(velocityField),
    g(dict.lookup("gravity")),
    T0(dict.lookup("T0")),
    nLayers(dict.lookup("BruntVaisalaFrequencies")),
    zLayers(dict.lookup("verticalLayers"))
{
    checkLayerSizes();
}

Foam::stratifiedThermalField::stratifiedThermalField
(
        const dimensionedVector& g,
        const dimensionedScalar& T0,
        const scalarList& bruntVaisalaFrequencies,
        const scalarList& verticalLayers,
        const advectable& velocityField
)
:
    tracerField(velocityField),
    g(g),
    T0(T0),
    nLayers(bruntVaisalaFrequencies),
    zLayers(verticalLayers)
{
    checkLayerSizes();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar Foam::stratifiedThermalField::tracerAt
(
        const Foam::point& p,
        const Foam::Time& t
) const
{
    const scalar z = p.z();

    scalar theta0 = T0.value();
    for(label il = 0; il < nLayers.size(); il++)
    {
        if (z >= zLayers[il] && z < zLayers[il+1])
        {
            return theta0*Foam::exp
            (
                sqr(nLayers[il])/mag(g.value())*(z-zLayers[il])
            );
        }
        theta0 *= Foam::exp
            (sqr(nLayers[il])/mag(g.value())*(zLayers[il+1]-zLayers[il]));
    }
    FatalErrorIn("stratifiedThermalField")
        << "height " << z << " not found in levels "
        << zLayers << " at point " << p << exit(FatalError);
    return -1;
}

vector Foam::stratifiedThermalField::gradAt
(
        const Foam::point& p,
        const Foam::Time& t
) const
{
    const scalar z = p.z();
    const scalar theta0 = T0.value();

    for (label il = 0; il < nLayers.size(); il++)
    {
        if (z >= zLayers[il] && z < zLayers[il+1])
        {
            const scalar N2 = sqr(nLayers[il]);
            const scalar gmag = mag(g).value();
            const scalar dtheta_dz = theta0 * N2 / gmag * Foam::exp(N2 / gmag * z);
            return vector(0, 0, dtheta_dz);
        }
    }
    FatalErrorIn("stratifiedThermalField")
        << "height " << z << " not found in levels "
        << zLayers << exit(FatalError);
    return vector(-1,-1,-1);
}

// ************************************************************************* //
