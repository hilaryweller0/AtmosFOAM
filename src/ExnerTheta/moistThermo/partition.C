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

#include "partition.H"
#include "moreListOps.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partition::partition
(
    const word& partitionName__,
    const wordList& partNames,
    const fvMesh& mesh,
    const dictionary dict,
    const volScalarField& Exner
)
:
    baseAtmosphere(add(partitionName__, partNames), mesh, dict),
    partitionName_(partitionName__),
    sigma_
    (
        IOobject
        (
            partitionName_+"sigma", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    sigmaRho_(sigma_*baseAtmosphere::sumDensity()),
    theta_
    (
        IOobject
        (
            partitionName_+"theta", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            partitionName_+"T", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        theta_*Exner
    ),
    Uf_
    (
        IOobject
        (
            partitionName_+"Uf", mesh.time().timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ),
    flux_
    (
        IOobject(partitionName_+"flux", mesh.time().timeName(), mesh),
        linearInterpolate(sigmaRho_)*(Uf_ & mesh.Sf())
    ),
    sigmaVolLiquid_
    (
        IOobject(partitionName_+"sigmaVolLiquid", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("sigmaVolLiquid", dimless, scalar(0))
    ),
    dSigmaRhodt_
    (
        IOobject(partitionName_+"dSigmaRhodt", mesh.time().timeName(), mesh),
        -fvc::div(flux_)
    ),
    dSigmaRhoThetadt_
    (
        IOobject(partitionName_+"dSigmaRhoThetadt", mesh.time().timeName(), mesh),
        -fvc::div(flux_, theta_)
    )
{
    sumVolLiquid();

    sigma_.oldTime();
    sigmaRho_.oldTime();
    dSigmaRhodt_.oldTime();
    theta_.oldTime();
    Uf_.oldTime();
    flux_.oldTime();
    dSigmaRhoThetadt_.oldTime();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partition::~partition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::partition::sumDensity()
{
    sigmaRho_ == dimensionedScalar("zero", dimDensity, scalar(0));
    for(label is = 0; is < size(); is++)
    {
        sigmaRho_ += operator[](is).gas().rho()
                  + operator[](is).liquid().rho()*operator[](is).liquid().v();
    }
    
    sigmaRho_ *= sigma_;
    
    return sigmaRho_;
}


Foam::volScalarField& Foam::partition::sumVolLiquid()
{
    sigmaVolLiquid_ = dimensionedScalar("zero", dimless, scalar(0));
    for(label is = 0; is < size(); is++)
    {
        sigmaVolLiquid_ += operator[](is).liquid().v();
    }
    sigmaVolLiquid_*= sigma_;
    return sigmaVolLiquid_;
}


Foam::volScalarField& Foam::partition::updateSigma(const volScalarField& p)
{
    sigma_ *= sumPressure(T_)/p;
    return sigma_;
}


//void Foam::partition::updateThetaT(const volScalarField& Exner)
//{
//    // Update theta and T from Exner using the equation of state for the
//    // perfect gas part of the partition
//    perfectGasPhase& air = operator[](0).gas();
//    const scalar kappa = air.kappa();
//    const dimensionedScalar& p0 = air.p0();
//    theta_ = p0*volGas()*pow(Exner, (1-kappa)/kappa)/rhoR();
//    T_ = theta_*Exner;
//}


void Foam::partition::write()
{
    baseAtmosphere::write();
    sigma_.write();
    sigmaRho_.write();
    T_.write();
    theta_.write();
    Uf_.write();
}


void Foam::partition::readUpdate(const volScalarField& Exner)
{
    const fvMesh& mesh = Exner.mesh();
    baseAtmosphere::readUpdate();
    sigma_ = volScalarField
    (
        IOobject(partitionName_+"sigma", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        sigma_
    );
    theta_ = volScalarField
    (
        IOobject(partitionName_+"theta", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        theta_
    );
    Uf_ = surfaceVectorField
    (
        IOobject(partitionName_+"Uf", mesh.time().timeName(), mesh,
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        Uf_
    );
    flux_ = linearInterpolate(sumDensity())*(Uf_ & mesh.Sf());
    sumDensity();
    sumVolLiquid();
    T_ = theta_*Exner;
}

// ************************************************************************* //
