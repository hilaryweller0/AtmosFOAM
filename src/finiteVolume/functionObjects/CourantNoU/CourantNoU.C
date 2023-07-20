/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "CourantNoU.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CourantNoU, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        CourantNoU,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::functionObjects::CourantNoU::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& U =
            lookupObject<volVectorField>(fieldName_);

        tmp<volScalarField> tCo
        (
            volScalarField::New
            (
                resultName_,
                mesh_,
                dimensionedScalar(dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        tCo->ref() = (0.5*mesh_.time().deltaT())/mesh_.V()*fvc::surfaceSum(mag
        (
            linearInterpolate(U) & mesh_.Sf()
        ))()();

        tCo->correctBoundaryConditions();

        return store(resultName_, tCo);
    }
    else
    {
        cannotFindObject<volVectorField>(fieldName_);

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNoU::CourantNoU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "Co", "U")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNoU::~CourantNoU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CourantNoU::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    return true;
}


// ************************************************************************* //
