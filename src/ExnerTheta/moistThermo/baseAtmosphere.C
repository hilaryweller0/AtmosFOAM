/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "baseAtmosphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::baseAtmosphere::baseAtmosphere
(
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict
)
:
    PtrList<fluidSpecie>(partNames.size())
{
    for(label ip = 0; ip < size(); ip++)
    {
        set
        (
            ip,
            new fluidSpecie
            (
                partNames[ip],
                IOobject(partNames[ip]+"VapourRho", mesh.time().timeName(), mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                IOobject(partNames[ip]+"LiquidFrac", mesh.time().timeName(), mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                mesh,
                dict.subDict(partNames[ip])
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::baseAtmosphere::~baseAtmosphere()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::volGas() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tvG
    (
        new volScalarField
        (
            IOobject
            (
                "volGas",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            1 - operator[](0).liquid().v()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        tvG.ref() -= operator[](ip).liquid().v();
    }
    return tvG;
}

Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::sumDensity() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tRho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
        )
    );

    for(label ip = 1; ip < size(); ip++)
    {
        tRho.ref() += operator[](ip).gas().rho()
                   + operator[](ip).liquid().rho()*operator[](ip).liquid().v();
    }
    return tRho;
}


Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::sumPressure
(
    const volScalarField& T
) const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> tp
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.partialPressure(T)
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        tp.ref() += operator[](ip).gas().partialPressure(T);
    }
    // Divide by the volume fraction occupied by gas
    tp.ref() /= volGas();
    return tp;
}


Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::rhoR() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoRt
    (
        new volScalarField
        (
            IOobject
            (
                "rhoRt",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.R()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoRt.ref() += operator[](ip).gas().rho()*operator[](ip).gas().R();
    }
    return trhoRt;
}

Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::rhoCp() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoCp
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCp",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.Cp()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
            *operator[](0).liquid().Cp()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoCp.ref() += operator[](ip).gas().rho()*operator[](ip).gas().Cp()
                    + operator[](ip).liquid().rho()*operator[](ip).liquid().v()
                      *operator[](ip).liquid().Cp();
    }
    return trhoCp;
}

Foam::tmp<Foam::volScalarField> Foam::baseAtmosphere::rhoCv() const
{
    const perfectGasPhase& air = operator[](0).gas();

    tmp<volScalarField> trhoCv
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCv",
                air.rho().time().timeName(),
                air.rho().mesh()
            ),
            air.rho()*air.Cv()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
            *operator[](0).liquid().Cp()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoCv.ref() += operator[](ip).gas().rho()*operator[](ip).gas().Cv()
                    + operator[](ip).liquid().rho()*operator[](ip).liquid().v()
                      *operator[](ip).liquid().Cp();
    }
    return trhoCv;
}


void Foam::baseAtmosphere::write()
{
    for(label ip = 0; ip < size(); ip++)
    {
        const fluidSpecie& phase = operator[](ip);
        phase.gas().rho().write();
        phase.liquid().v().write();
    }
}


void Foam::baseAtmosphere::readUpdate()
{
    for(label ip = 0; ip < size(); ip++)
    {
        fluidSpecie& phase = operator[](ip);
        const fvMesh& mesh = phase.gas().rho().mesh();
        
        phase.gas().rho() = volScalarField
        (
            IOobject(phase.name()+"VapourRho", mesh.time().timeName(), mesh,
                     IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            phase.gas().rho()
        );
        phase.liquid().v() = volScalarField
        (
            IOobject(phase.name()+"LiquidFrac", mesh.time().timeName(), mesh,
                     IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            phase.liquid().v()
        );
    }
}


// ************************************************************************* //
