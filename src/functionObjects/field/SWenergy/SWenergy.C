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

#include "SWenergy.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(SWenergy, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        SWenergy,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::SWenergy::writeFileHeader(const label i)
{
    writeHeader(file(), "SWenergy ()");
    writeCommented(file(), "Time");
    writeTabbed(file(), "average");
    writeTabbed(file(), "max");
    file() << endl;
}


void Foam::functionObjects::SWenergy::calcSWenergy()
{
    const volScalarField& h = mesh_.lookupObject<volScalarField>("h");
    const volScalarField& h0 = mesh_.lookupObject<volScalarField>("h0");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const IOdictionary& earthP
         = mesh_.lookupObject<IOdictionary>("earthProperties");
    const dimensionedScalar magg(earthP.lookup("magg"));
    
    KE_ = 0.5*h*magSqr(U);
    PE_ = magg*h*(0.5*h + h0);
    SWenergy_ = KE_ + PE_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::SWenergy::SWenergy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    KE_
    (
        IOobject("KE", mesh_),
        mesh_,
        dimensionedScalar(dimensionSet(0,3,-2,0,0), 0)
    ),
    PE_
    (
        IOobject("PE", mesh_),
        mesh_,
        dimensionedScalar(dimensionSet(0,3,-2,0,0), 0)
    ),
    SWenergy_
    (
        IOobject("SWenergy", mesh_),
        mesh_,
        dimensionedScalar(dimensionSet(0,3,-2,0,0), 0)
    )
{
    calcSWenergy();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::SWenergy::~SWenergy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::SWenergy::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    return true;
}


bool Foam::functionObjects::SWenergy::execute()
{
    calcSWenergy();
    logFiles::write();

    if (Pstream::master())
    {
        scalar energyMax = gMax(SWenergy_.primitiveField());
        scalar energyMean
             = gSum(SWenergy_.primitiveField()*mesh_.V().primitiveField())
                    /gSum(mesh_.V().primitiveField());

        scalar KEMax = gMax(KE_.primitiveField());
        scalar KEMean
             = gSum(KE_.primitiveField()*mesh_.V().primitiveField())
                    /gSum(mesh_.V().primitiveField());
    
        scalar PEMax = gMax(PE_.primitiveField());
        scalar PEMean
             = gSum(PE_.primitiveField()*mesh_.V().primitiveField())
                    /gSum(mesh_.V().primitiveField());
    
        Log << "Time "
            << mesh_.time().timeName(mesh_.time().userTimeValue())
            << " Energy max: " << energyMax << " average: " << energyMean
            << " KE max: " << KEMax << " average: " << KEMean
            << " PE max: " << PEMax << " average: " << PEMean
             << nl;

/*        writeTime(file());
        file()
            << tab << energyMean
            << tab << energyMax
            << endl;*/
    }

    return true;
}


bool Foam::functionObjects::SWenergy::write()
{
    writeLocalObjects::write();
    SWenergy_.write();
    KE_.write();
    PE_.write();
    logFiles::write();
    return true;
}


// ************************************************************************* //
