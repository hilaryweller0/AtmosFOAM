/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2025 OpenFOAM Foundation
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

#include "theta.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(theta, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        theta,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::theta::writeFileHeader(const label i)
{
    writeHeader(file(), "theta ()");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


Foam::volScalarField& Foam::functionObjects::theta::calcTheta()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    theta_ = T*pow(pRef_/p, RbyCp());
    return theta_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::theta::theta
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    Tref_(dimTemperature, readScalar(dict.lookup("Tref"))),
    pRef_(mesh_.lookupObject<uniformDimensionedScalarField>("pRef")),
    RbyCp_(0),
    theta_
    (
        IOobject(IOobject::groupName(type(), phaseName_), mesh_),
        mesh_,
        Tref_
    )
{
    IOdictionary
        thermoDict(mesh_.lookupObject<IOdictionary>("physicalProperties"));
    scalar W = readScalar
    (
        thermoDict.subDict("mixture").subDict("specie").lookup("molWeight")
    );
    scalar Cp = readScalar
    (
        thermoDict.subDict("mixture").subDict("thermodynamics").lookup("Cp")
    );
    RbyCp_ = constant::thermodynamic::RR/(W*Cp);
    calcTheta();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::theta::~theta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::theta::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    resetName(IOobject::groupName(typeName, phaseName_));
    resetLocalObjectName(IOobject::groupName(typeName, phaseName_));
    return true;
}


bool Foam::functionObjects::theta::execute()
{
    calcTheta();
    logFiles::write();
    return true;
}


bool Foam::functionObjects::theta::write()
{
    //Log << type() << " " << name() << " write:" << nl;
    writeLocalObjects::write();
    theta_.write();
    logFiles::write();
    return true;
}


// ************************************************************************* //
