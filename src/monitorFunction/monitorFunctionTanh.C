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

#include "monitorFunctionTanh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(monitorFunctionTanh, 0);
addToRunTimeSelectionTable(monitorFunction, monitorFunctionTanh, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

monitorFunctionTanh::monitorFunctionTanh
(
    const IOdictionary& dict
)
:
    monitorFunction(dict),
    alpha_(readScalar(dict.lookup("alpha"))),
    beta_(readScalar(dict.lookup("beta"))),
    gamma_(readScalar(dict.lookup("gamma"))),
    centre_(dict.lookup("centre")),
    spherical_(dict.lookup("spherical"))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

//tmp<scalarField> monitorFunctionTanh::map
//(
//    const pointField& pts,
//    const scalarField& oldMonitor
//) const
//{
//    tmp<scalarField> tMon(new scalarField(pts.size(), scalar(0)));
//    scalarField& mon = tMon();

//    point centreHat = spherical_? point(unitVector(centre_)) : centre_;
//    
//    forAll(mon, ip)
//    {
//        scalar r = spherical_? acos
//                  (
//                    point(unitVector(pts[ip])) & centreHat
//                  )
//                 : mag(pts[ip] - centre_);
//        mon[ip] = sqrt(0.5/(1+gamma_)*(tanh((beta_-r)/alpha_)+1) + gamma_);
//    }
//    
//    return tMon;
//}


tmp<volScalarField> monitorFunctionTanh::map
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
            dimensionSet(0,0,0,0,0),
            wordList(newMesh.boundaryMesh().size(), "zeroGradient")
        )
    );
    volScalarField& mon = tMon();
    
    point centreHat = spherical_? point(unitVector(centre_)) : centre_;
    
    forAll(mon, cellI)
    {
        scalar r = spherical_? acos
                  (
                    point(unitVector(newMesh.C()[cellI])) & centreHat
                  )
                 : mag(newMesh.C()[cellI] - centre_);
        mon[cellI] = sqrt(0.5/(1+gamma_)*(tanh((beta_-r)/alpha_)+1) + gamma_);
    }
    
    mon.correctBoundaryConditions();
    return tMon;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
