/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "uStar.H"
#include "momentumTransportModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(uStar, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        uStar,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::uStar::writeFileHeader(const label i)
{
    writeHeader(file(), "uStar ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::uStar::calcuStar
(
    const momentumTransportModel& turbModel
)
{
    tmp<volScalarField> tuStar
    (
        volScalarField::New
        (
            IOobject::groupName(type(), phaseName_),
            mesh_,
            dimensionedScalar(dimVelocity, 0)
        )
    );

    volScalarField::Boundary& uStarBf = tuStar.ref().boundaryFieldRef();

    const tmp<surfaceVectorField> tdevTau = turbModel.devTau();
    const surfaceVectorField::Boundary& devTauBf = tdevTau().boundaryField();

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            uStarBf[patchi] = sqrt(mag(devTauBf[patchi]));
        }
    }

    return tuStar;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::uStar::uStar
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::uStar::~uStar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::uStar::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    resetName(IOobject::groupName(typeName, phaseName_));
    resetLocalObjectName(IOobject::groupName(typeName, phaseName_));

    return true;
}


bool Foam::functionObjects::uStar::execute()
{
    if (mesh_.foundObject<momentumTransportModel>
    (
        IOobject::groupName(momentumTransportModel::typeName, phaseName_))
    )
    {
        const momentumTransportModel& model =
            mesh_.lookupType<momentumTransportModel>(phaseName_);

        store(IOobject::groupName(type(), phaseName_), calcuStar(model));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    const volScalarField& uStar =
        mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(type(), phaseName_)
        );

    const volScalarField::Boundary& uStarBf = uStar.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    logFiles::write();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& uStarp = uStarBf[patchi];

            const scalar minuStar = gMin(uStarp);
            const scalar maxuStar = gMax(uStarp);
            const scalar avguStar = gAverage(uStarp);

            if (Pstream::master())
            {
                Log << "Time "
                    << mesh_.time().timeName(mesh_.time().userTimeValue())
                    << " patch " << patch.name()
                    << " uStar: min: " << minuStar << " max: " << maxuStar
                    << " average: " << avguStar << nl;

                writeTime(file());
                file()
                    << tab << patch.name()
                    << tab << minuStar
                    << tab << maxuStar
                    << tab << avguStar
                    << endl;
            }
        }
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::uStar::write()
{
    Log << type() << " " << name() << " write:" << nl;
    writeLocalObjects::write();
    logFiles::write();
    return true;
}


// ************************************************************************* //
