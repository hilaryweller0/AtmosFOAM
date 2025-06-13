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

#include "multiScalarRKTransport.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcFlux.H"
#include "fvModels.H"

#include "EulerDdtScheme.H"
#include "localMax.H"
#include "zeroGradientFvPatchFields.H"
#include "CourantNoFunc.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(multiScalarRKTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        multiScalarRKTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::multiScalarRKTransport::multiScalarRKTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    nFields_(readLabel(dict.lookup("nFields"))),
    fieldNames_(dict.lookup("fields")),
    s_(nFields_),
    phivName_(dict.lookup("advectingFlux")),
    massFluxName_(dict.lookup("massFlux")),
    RK_(dict.lookup("RK_ButcherCoeffs")),
    gammaCoeffs_(dict.lookup("gammaCoeffs")),
    correctionSchemeName_(dict.lookup("correctionScheme"))
{
    if (nFields_ != fieldNames_.size())
    {
        FatalErrorIn("multiScalarRKTransport::multiScalarRKTransport")
            << "nFields = " << nFields_ << " but fields are " << fieldNames_
            << exit(FatalError);
    }
    
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::multiScalarRKTransport::~multiScalarRKTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::multiScalarRKTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    for(label is = 0; is < nFields_; is++)
    {
        Info << "Reading " << fieldNames_[is] << endl;
        s_.set
        (
            is,
            volScalarField
            (
                IOobject
                (
                    fieldNames_[is], time_.name(), mesh_, IOobject::MUST_READ
                ),
                mesh_
            )
        );

        //s_[is].writeOpt() = IOobject::AUTO_WRITE;
        s_[is].checkIn();
        s_[is].rename(fieldNames_[is]);
        s_[is].correctBoundaryConditions();
        s_[is].oldTime();
    }
    
    return true;
}


Foam::wordList Foam::functionObjects::multiScalarRKTransport::fields() const
{
    return fieldNames_;
}


bool Foam::functionObjects::multiScalarRKTransport::execute()
{
    Info<< type() << " execute:" << endl;

    // Reference to the time step
    const dimensionedScalar& dt = time_.deltaT();

    // Look up or read the advectingFlux, phiv
    tmp<surfaceScalarField> phivt;
    if(mesh_.foundObject<surfaceScalarField>(phivName_))
    {
        Info << "Looking up " << phivName_ << endl;
        phivt = &mesh_.lookupObjectRef<surfaceScalarField>(phivName_);
        Info << "Done" << endl;
    }
    else
    {
        Info << "Reading " << phivName_ << endl;
        phivt = new surfaceScalarField
        (
            IOobject
            (
                phivName_,
                s_[0].time().startTime().name(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
    }
    surfaceScalarField& phiv = phivt.ref();
    
    // Initialise the necessary schemes
    localMax<scalar> maxInterp(mesh_);
    fv::EulerDdtScheme<scalar> backwardEuler(mesh_);
    fv::gaussConvectionScheme<scalar> upwindCon0
        (mesh_, phiv.oldTime(), upwindScheme(phiv.oldTime()));
    fv::gaussConvectionScheme<scalar> upwindCon(mesh_, phiv, upwindScheme(phiv));

    // Calculate the maximum Courant number on faces
    volScalarField Co0 = CourantNo(phiv.oldTime(), dt);
    volScalarField Co1 = CourantNo(phiv, dt);
    surfaceScalarField Co = max
    (
        maxInterp.interpolate(Co0), maxInterp.interpolate(Co1)
    );
    // Calculate alpha, beta and gamma as a function of the max Courant number
    surfaceScalarField alpha = 1 - 1/max(scalar(2), Co);
    surfaceScalarField beta  = 1 - 1/max(scalar(1), Co);
    surfaceScalarField gamma = gammaCoeffs_[1]
                /(gammaCoeffs_[1]-gammaCoeffs_[0] + max(Co, gammaCoeffs_[0]));
    Info << "Co goes from " << min(Co).value() << " to " << max(Co).value() << nl
         << "alpha goes from " << min(alpha).value() << " to " << max(alpha).value() << nl
         << "beta goes from " << min(beta).value() << " to " << max(beta).value() << nl
         << "gamma goes from " << min(gamma).value() << " to " << max(gamma).value() 
         << endl;

    // Accumulate the total flux
    PtrList<surfaceScalarField> totalFlux(nFields_);
    
    // Total fluxes for the RK stages for all tracers
    PtrList<PtrList<surfaceScalarField>> F(RK_.n()+1);
    for(label iRK = 0; iRK <= RK_.n(); iRK++)
    {
        F.set(iRK, new PtrList<surfaceScalarField>(nFields_));
    }
    
    // The density for each RK stage
    PtrList<volScalarField> rho(RK_.n()+1);
    if (massFluxName_ != "none") rho.set(0, s_[0]);
    
    // The zeroth flux is for the explicit low order update
    surfaceScalarField phi = phiv.oldTime();
    for(label is = 0; is < nFields_; is++)
    {
        // Is this field density weighted?
        const bool densityWeighted = is > 0 && massFluxName_ != "none";
        
        // Calculate the first low-order flux
        F[0].set
        (
            is,
            (1-alpha)*beta*phi*upInterp(phi,s_[is])
        );
        // Accumulate the total flux, starting from the old time step terms
        totalFlux.set(is, F[0][is]);
        
        // Advance s by the low-order flux
        Info << "Updating " << s_[is].name() << " by maximum "
                << max(mag(fvc::div(F[0][is]))) << endl;
        if (densityWeighted)
        {
            s_[is] = (s_[is].oldTime()*rho[0] - dt*fvc::div(F[0][is]))/rho[1];
        }
        else
        {
            s_[is] -= dt*fvc::div(F[0][is]);
        }
        
        // Change the flux for following, non-density tracers
        if (is == 0 && massFluxName_ != "none")
        {
            phi *= upInterp(phi,s_[0]);
            rho.set(1, s_[0]);
        }
    }
    
    // RK advection
    for(int iRK = 0; iRK < RK_.n(); iRK++)
    {
        // Sub-stage size
        scalar c = 0;
        for(int j = 0; j <= iRK; j++) c += RK_[iRK][j];
        
        // Advecting flux at the sub time
        phi = (1-c)*phiv.oldTime() + c*phiv;

    }
    
    return true;
}


bool Foam::functionObjects::multiScalarRKTransport::write()
{
    for(label is = 0; is < nFields_; is++)
    {
        s_[is].write();
    }
    return true;
}


// ************************************************************************* //
