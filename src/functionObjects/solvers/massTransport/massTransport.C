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

#include "massTransport.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcFlux.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"

#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "localMax.H"
#include "zeroGradientFvPatchFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        massTransport,
        dictionary
    );
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::massTransport::diffusivityType,
    3
>::names[] =
{
    "none",
    "constant",
    "viscosity"
};

const Foam::NamedEnum
<
    Foam::functionObjects::massTransport::diffusivityType,
    3
> Foam::functionObjects::massTransport::diffusivityTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::massTransport::D() const
{
    const word Dname("D" + fieldName_);

    if (diffusivity_ == diffusivityType::constant)
    {
        return volScalarField::New
        (
            Dname,
            mesh_,
            dimensionedScalar(Dname, dimKinematicViscosity, D_)
        );
    }
    else
    {
        const momentumTransportModel& turbulence =
            mesh_.lookupType<momentumTransportModel>();

        return volScalarField::New
        (
            Dname,
            alphal_*turbulence.nu() + alphat_*turbulence.nut()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massTransport::massTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")),
    diffusivity_(diffusivityType::none),
    D_(0),
    nCorr_(0),
    s_
    (
        IOobject
        (
            fieldName_,
            time_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massTransport::~massTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("advectingFlux", "phi");
    fluxName_ = dict.lookupOrDefault<word>("scalarFlux", "none");
    rhoName_ = dict.lookupOrDefault<word>("density", "none");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_);
    diffusivity_ = diffusivityTypeNames_.read(dict.lookup("diffusivity"));

    switch(diffusivity_)
    {
        case diffusivityType::none:
            break;

        case diffusivityType::constant:
            dict.lookup("D") >> D_;
            break;

        case diffusivityType::viscosity:
            dict.lookup("alphal") >> alphal_;
            dict.lookup("alphat") >> alphat_;
            break;
    }

    nCorr_ = dict.lookupOrDefaultBackwardsCompatible<label>
    (
        {"nCorrectors", "nCorr"},
        0
    );

    return true;
}


Foam::wordList Foam::functionObjects::massTransport::fields() const
{
    return wordList{phiName_};
}


bool Foam::functionObjects::massTransport::execute()
{
    Info<< type() << " execute:" << endl;

    // Look up or read the advectingFlux, phiv
    tmp<surfaceScalarField> phit;
    if(mesh_.foundObject<surfaceScalarField>(phiName_))
    {
        Info << "Looking up " << phiName_ << endl;
        phit = &mesh_.lookupObjectRef<surfaceScalarField>(phiName_);
        Info << "Done" << endl;
    }
    else
    {
        Info << "Reading in " << phiName_ << endl;
        phit = new surfaceScalarField
        (
            IOobject
            (
                phiName_,
                s_.time().startTime().name(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );
        Info << "done" << endl;
    }
    surfaceScalarField& phi = phit.ref();
    
    const word divScheme("div(" + phiName_ + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.solution().relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.solution().equationRelaxationFactor(schemesField_);
    }

    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(mesh_)
    );

    if (phi.dimensions() == dimVolume/dimTime)
    {
        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s_)
              + fvm::div(phi, s_, divScheme)
             ==
                fvModels.source(s_)
            );

            if (diffusivity_ != diffusivityType::none)
            {
                sEqn -= fvm::laplacian(D(), s_);
            }

            sEqn.relax(relaxCoeff);

            fvConstraints.constrain(sEqn);

            sEqn.solve(schemesField_);

            fvConstraints.constrain(s_);
            
            if (fluxName_ != "none") // Update the scalarFlux in the database
            {
                Info << "Recalculating " << fluxName_ << " after advection of " 
                     << fieldName_ << endl;
                surfaceScalarField& flux =
                    mesh_.lookupObjectRef<surfaceScalarField>(fluxName_);
                flux = sEqn.flux();
            }
        }
        const surfaceScalarField& fluxTmp 
            = mesh_.lookupObject<surfaceScalarField>(fluxName_);
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s_)
              + fvm::div(phi, s_, divScheme)
             ==
                fvModels.source(rho, s_)
            );

            if (diffusivity_ != diffusivityType::none)
            {
                sEqn -= fvm::laplacian(rho*D(), s_);
            }

            sEqn.relax(relaxCoeff);

            fvConstraints.constrain(sEqn);

            sEqn.solve(schemesField_);

            fvConstraints.constrain(s_);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::massTransport::write()
{
    s_.write();
    return true;
}


// ************************************************************************* //
