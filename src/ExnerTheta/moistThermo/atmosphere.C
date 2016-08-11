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

#include "atmosphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphere::atmosphere
(
    const wordList& partNames,
    const word timeName,
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
                IOobject(partNames[ip]+"VapourRho", timeName, mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                IOobject(partNames[ip]+"LiquidFrac", timeName, mesh,
                         IOobject::MUST_READ, IOobject::AUTO_WRITE),
                mesh,
                dict.subDict(partNames[ip])
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphere::~atmosphere()
{
 The constructor calls new so I think that the destructor should do something
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::atmosphere::sumDensity()
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rho", 
                operator[](0).gas().rho().time().timeName(),
                operator[](0).gas().rho().mesh()
            ),
            operator[](0).gas().rho()
          + operator[](0).liquid().rho()*operator[](0).liquid().v()
        )
    );
    for(label ip = 1; ip < size(); ip++)
    {
        trhoRef() += operator[](ip).gas().rho()
                  + operator[](ip).liquid().rho()*operator[](0).liquid().v();
    }
    return trho;
    
}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::gasConst()
{

}

Foam::tmp<Foam::volScalarField> Foam::atmosphere::Cp()

Foam::tmp<Foam::volScalarField> Foam::atmosphere::Cv()

Foam::tmp<Foam::volScalarField> Foam::atmosphere::sumPressure()

Foam::tmp<Foam::volScalarField> Foam::atmosphere::thetaRho
(
    const volScalarField& theta
)

Foam::tmp<Foam::volScalarField> Foam::atmosphere::pFromT(const volScalarField& T)

Foam::tmp<Foam::volScalarField> Foam::atmosphere::ExnerFromTheta
(
    const volScalarField& theta
)

Foam::tmp<Foam::volScalarField> Foam::atmosphere::ExnerFromThetaRho
(
    const volScalarField& thetaRho
)

Foam::tmp<Foam::volScalarField> Foam::atmosphere::thetaSource
(
    const volScalarField& T, const volScalarField& Scond
)


// ************************************************************************* //
