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

#include "axialAngMom.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(axialAngMom, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        axialAngMom,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::axialAngMom::writeFileHeader(const label i)
{
    writeHeader(file(), "axialAngMom ()");
    writeCommented(file(), "Time");
    writeTabbed(file(), "average");
    writeTabbed(file(), "normalised");
    file() << endl;
}


void Foam::functionObjects::axialAngMom::calcAxialAngMom()
{
    const volScalarField& h = mesh_.lookupObject<volScalarField>("h");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    
    axialAngMom_ = h*(U & east_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::axialAngMom::axialAngMom
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    east_(vector(0,0,1) ^ mesh_.C()),
    axialAngMom0_(dimensionSet(0,3,-1,0,0), 0),
    axialAngMom_
    (
        IOobject("axialAngMom", mesh_),
        mesh_,
        axialAngMom0_
    )
{
    calcAxialAngMom();
    axialAngMom0_ = sum(axialAngMom_*mesh_.V())
                    /sum(mesh_.V());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::axialAngMom::~axialAngMom()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::axialAngMom::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    return true;
}


bool Foam::functionObjects::axialAngMom::execute()
{
    calcAxialAngMom();
    logFiles::write();

    if (Pstream::master())
    {
        scalar momMean
             = gSum(axialAngMom_.primitiveField()*mesh_.V().primitiveField())
                    /gSum(mesh_.V().primitiveField());
    
        Log << "Axial Angular Momentum at time "
            << mesh_.time().timeName(mesh_.time().userTimeValue())
            << " average: " << momMean
            << " normalised: "
            << (momMean - axialAngMom0_.value())/axialAngMom0_.value()
             << nl;

/*        writeTime(file());
        file()
            << tab << momMean
            << tab << (momMean - axialAngMom0_)/axialAngMom0_
            << endl;*/
    }

    return true;
}


bool Foam::functionObjects::axialAngMom::write()
{
    writeLocalObjects::write();
    axialAngMom_.write();
    logFiles::write();
    return true;
}


// ************************************************************************* //
