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

#include "heatFlux.H"
#include "rhoThermo.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "ThermophysicalTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "physicalProperties.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "constTransport.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(heatFlux, 0);
    addToRunTimeSelectionTable(functionObject, heatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::heatFlux::calc()
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
        const volScalarField& T = lookupObject<volScalarField>(TName_);

        Info<< "Creating turbulence model\n" << endl;
        autoPtr<compressible::momentumTransportModel> turbulence
        (
            compressible::momentumTransportModel::New(rho, U, phi, thermo)
        );

        Info<< "Creating thermophysical transport model\n" << endl;
        typedef ThermophysicalTransportModel
        <
            compressible::momentumTransportModel, fluidThermo
        > tModel;

        autoPtr<tModel> thermophysicalTransport(tModel::New(turbulence(), thermo));

        const dimensionedScalar Tref
        (
            "Tref", dimTemperature, readScalar(thermo.properties().lookup("Tref"))
        );
        const dimensionedScalar pRef
        (
            "pRef", dimPressure, readScalar(thermo.properties().lookup("pRef"))
        );
        const constTransport<hConstThermo<perfectGas<specie> > > air
        (
            "mixture", thermo.properties()
        );

        const dimensionedScalar Cp
        (
            "Cp", dimGasConstant, air.Cp(pRef.value(),Tref.value())
        );
        
        surfaceScalarField alphaf
             = linearInterpolate(thermophysicalTransport->kappaEff()/Cp);
        surfaceScalarField Tf = linearInterpolate(T);
        surfaceScalarField qFlux = phi*Tf - alphaf*mesh_.magSf()*fvc::snGrad(T);

        tmp<volVectorField> tq
        (
            new volVectorField
            (
                IOobject(resultName_, mesh_.time().timeName(), mesh_),
                mesh_,
                dimensionedVector
                (
                    rho.dimensions()*U.dimensions()*T.dimensions(),
                    vector::zero
                ),
                "calculated"
            )
        );
        volVectorField& q = tq.ref();
        q = fvc::reconstruct(qFlux);
        
        forAll(q.boundaryFieldRef(), patchi)
        {
            q.boundaryFieldRef()[patchi] = qFlux.boundaryField()[patchi]
                    *mesh_.Sf().boundaryField()[patchi]
                         / sqr(mesh_.magSf().boundaryField()[patchi]);
        }
        
        return store(resultName_, tq);
    }
    else
    {
        cannotFindObject<surfaceScalarField>(fluxName_);

        return false;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::heatFlux::heatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldsExpression(name, runTime, dict)
{
    read(dict);
    resultName_ = "heatFlux";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::heatFlux::~heatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::heatFlux::read(const dictionary& dict)
{
    dict.readIfPresent("U", UName_);
    dict.readIfPresent("rho", rhoName_);
    dict.readIfPresent("phi", fluxName_);
    dict.readIfPresent("T", TName_);
    return true;
}

// ************************************************************************* //
