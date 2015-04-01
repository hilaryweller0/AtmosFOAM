/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "monitorFunctionSech.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(monitorFunctionSech, 0);
addToRunTimeSelectionTable(monitorFunction, monitorFunctionSech, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

monitorFunctionSech::monitorFunctionSech
(
    const IOdictionary& dict
)
:
    monitorFunction(dict),
    alpha1_(readScalar(dict.lookup("alpha1"))),
    alpha2_(readScalar(dict.lookup("alpha2"))),
    a_(readScalar(dict.lookup("a"))),
    sqra_(sqr(a_)),
    centre_(dict.lookup("centre"))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

//tmp<scalarField> monitorFunctionSech::map
//(
//    const pointField& pts,
//    const scalarField& oldMonitor
//) const
//{
//    tmp<scalarField> tMon(new scalarField(pts.size(), scalar(0)));
//    scalarField& mon = tMon();

//    forAll(mon, ip)
//    {
//        scalar R = mag(pts[ip] - centre_);
//        mon[ip] = 1 + alpha1_/sqr(cosh(alpha2_*(sqr(R) - sqra_)));
//    }
//    
//    return tMon;
//}


tmp<volScalarField> monitorFunctionSech::map
(
    const fvMesh& newMesh,
    const volScalarField& oldMonitor
) const
{
    tmp<volScalarField> tMon
    (
        new volScalarField
        (
            IOobject("monitor", newMesh.time().timeName(), newMesh,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            newMesh,
            dimensionSet(0,-2,0,0,0),
            wordList(newMesh.boundaryMesh().size(), "zeroGradient")
        )
    );
    volScalarField& mon = tMon();
    
    forAll(mon, cellI)
    {
        scalar R = mag(newMesh.C()[cellI] - centre_);
        mon[cellI] = 1 + alpha1_/sqr(cosh(alpha2_*(sqr(R) - sqra_)));
    }
    
    mon.correctBoundaryConditions();
    return tMon;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
