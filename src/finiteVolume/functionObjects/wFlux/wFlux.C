/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "wFlux.H"
#include "rhoThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wFlux, 0);
    addToRunTimeSelectionTable(functionObject, wFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::wFlux::calc()
{
    // The thermo model
    autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh_));
    rhoThermo& thermo = pThermo();

    // Read in fields
    if (foundObject<surfaceScalarField>(fluxName_))
    {
        const volVectorField& U = lookupObject<volVectorField>(UName_);
        const surfaceScalarField& phi = lookupObject<surfaceScalarField>(fluxName_);
        const volScalarField& rho = lookupObject<volScalarField>(rhoName_);
        const volScalarField& p = lookupObject<volScalarField>(pName_);

        const dimensionedScalar pRef = fvc::domainIntegrate(p)
                /fvc::domainIntegrate(p/p);
        //(
        //  "pRef", dimPressure, readScalar(thermo.properties().lookup("pRef"))
        //);

        autoPtr<compressible::momentumTransportModel> turbulence
        (
            compressible::momentumTransportModel::New(rho, U, phi, thermo)
        );

        surfaceScalarField muf = linearInterpolate(rho*turbulence->nuEff());
        surfaceScalarField wf = linearInterpolate(U.component(2));
        volVectorField wFlux = fvc::reconstruct
        (
            phi*wf - mesh_.magSf()*muf*fvc::snGrad(U.component(2))
        );
        tmp<volScalarField> tFlux(wFlux.component(2)+ (p - pRef));

        return store(resultName_, tFlux);
    }
    else
    {
        cannotFindObject<surfaceScalarField>(fluxName_);

        return false;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wFlux::wFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldsExpression(name, runTime, dict)
{
    read(dict);
    resultName_ = "wFlux";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wFlux::~wFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wFlux::read(const dictionary& dict)
{
    dict.readIfPresent("U", UName_);
    dict.readIfPresent("rho", rhoName_);
    dict.readIfPresent("phi", fluxName_);
    dict.readIfPresent("p", pName_);
    return true;
}

// ************************************************************************* //
