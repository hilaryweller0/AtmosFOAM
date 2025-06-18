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
#include "fvcFluxLimit.H"

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
    correctionSchemeName_(dict.lookup("correctionScheme")),
    FCTiter_(dict.lookup("FCTiters"))
{
    if (nFields_ != fieldNames_.size() || nFields_ != FCTiter_.size())
    {
        FatalErrorIn("multiScalarRKTransport::multiScalarRKTransport")
            << "nFields = " << nFields_ << " but fields are " << fieldNames_
            << " and FCTiters are " << FCTiter_
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
    const surfaceScalarField phiv 
      = mesh_.foundObject<surfaceScalarField>(phivName_) ?
        mesh_.lookupObjectRef<surfaceScalarField>(phivName_) :
        surfaceScalarField
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
    
    // Initialise the local max interpolation to find the Courant number
    localMax<scalar> maxInterp(mesh_);

    // Calculate the maximum Courant number on faces
    surfaceScalarField Co = max
    (
        maxInterp.interpolate(CourantNo(phiv.oldTime(), dt)),
        maxInterp.interpolate(CourantNo(phiv, dt))
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

    // Total fluxes for the RK stages for all tracers
    PtrList<PtrList<surfaceScalarField>> F(RK_.n()+1);
    for(label iRK = 0; iRK <= RK_.n(); iRK++)
    {
        F.set(iRK, new PtrList<surfaceScalarField>(nFields_));
    }
    
    // First Strang low-order update, before the RK stages, for all fields
    // First, the advecting flux is not density weighted
    surfaceScalarField phi = (1-alpha)*beta*phiv.oldTime();
    for(label is = 0; is < nFields_; is++)
    {
        // Calculate the first low-order flux
        F[0].set(is, phi*upInterp(phi,s_[is]));

        // Change the flux for following, non-density tracers
        if (is == 0 && massFluxName_ != "none")
        {
            phi *= upInterp(phi,s_[0]);
        }
    }
    for(label is = 0; is < nFields_; is++)
    {
        // Is this field density weighted?
        const bool densityWeighted = is > 0 && massFluxName_ != "none";
        // Advance s by the low-order flux
        if (!densityWeighted)
        {
            s_[is] -= dt*fvc::div(F[0][is]);
        }
        else
        {
            s_[is] = (s_[is].oldTime()*s_[0].oldTime() - dt*fvc::div(F[0][is]))
                     /s_[0];
        }
    }
    
    // Explicit RK advection stages
    for(int iRK = 0; iRK < RK_.n(); iRK++)
    {
        // Sub-stage size
        scalar c = 0;
        for(int j = 0; j <= iRK; j++) c += RK_[iRK][j];
        
        // Advecting flux at the sub time
        phi.dimensions().reset(phiv.dimensions());
        phi = (1-c)*phiv.oldTime() + c*phiv;

        // The flux for eac RK stage for each trancer
        for(int is = 0; is < nFields_; is++)
        {
            // Is this field density weighted?
            const bool densityWeighted = is > 0 && massFluxName_ != "none";
        
            // Calculate the scalar flux for this stage
            if (!densityWeighted)
            {
                F[iRK+1].set
                (
                    is,
                    phi*
                    (
                        (1-beta)*upInterp(phi,s_[is])
                      + gamma*hCorr(phi, s_[is])
                    )
                );
            }
            else
            {
                F[iRK+1].set
                (
                    is,
                    phi*
                    (
                        (1-beta)*upInterp(phi,s_[0]*s_[is])
                      + gamma*hCorr(phi, s_[0]*s_[is])
                    )
                );
            }
        }
        // This RK stage for each trancer
        for(int is = 0; is < nFields_; is++)
        {
            // Is this field density weighted?
            const bool densityWeighted = is > 0 && massFluxName_ != "none";

            //Update the scalar for each stage
            if (!densityWeighted)
            {
                s_[is] = s_[is].oldTime() - dt*fvc::div(F[0][is]);
                for(int j = 0; j <= iRK; j++)
                {
                    s_[is] -= dt*RK_[iRK][j]*fvc::div(F[j+1][is]);
                }
            }
            else
            {
                volScalarField rhos = s_[is].oldTime()*s_[0].oldTime()
                                    - dt*fvc::div(F[0][is]);
                for(int j = 0; j <= iRK; j++)
                {
                    rhos -= dt*RK_[iRK][j]*fvc::div(F[j+1][is]);
                }
                s_[is] = rhos/s_[0];
            }
        }
    }

    // Accumulate the total flux for each scalar
    PtrList<surfaceScalarField> totalFlux(nFields_);
    for(int is = 0; is < nFields_; is++)
    {
        totalFlux.set(is, F[0][is]);
        for(label j = 0; j < RK_.n(); j++)
        {
            totalFlux[is] += RK_[RK_.n()-1][j]*F[j+1][is];
        }
    }

    // Final, implicit, low-order Strang stage
    phi.dimensions().reset(phiv.dimensions());
    phi = alpha*beta*phiv;
    fv::EulerDdtScheme<scalar> EulerDdt(mesh_);
    fv::gaussConvectionScheme<scalar> upwindCon(mesh_, phi, upwindScheme(phi));
    for(label is = 0; is < nFields_; is++)
    {
        // Is this field density weighted?
        const bool densityWeighted = is > 0 && massFluxName_ != "none";
        
        // Create the matrix equation
        fvScalarMatrix sEqn
        (
            upwindCon.fvmDiv(phi, s_[is])
          + fvc::div(totalFlux[is])
        );
        
        // Add the rate of change term
        if (!densityWeighted)
        {
            sEqn += fvScalarMatrix(EulerDdt.fvmDdt(s_[is]));
        }
        else
        {
            sEqn += fvScalarMatrix(EulerDdt.fvmDdt(s_[0], s_[is]));
        }
        
        // Solve and update the total flux
        sEqn.solve();
        totalFlux[is] += sEqn.flux();
        
        // Change the flux for following, non-density tracers
        if (is == 0 && massFluxName_ != "none")
        {
            phi *= upInterp(phi,s_[0]);
        }

        // Apply FCT if needed
        if (FCTiter_[is] > 0)
        {
            surfaceScalarField fluxCorr = totalFlux[is];
            
            // Find the bounded solution and the bounded flux
            if (!densityWeighted)
            {
                const surfaceScalarField phi0 = (1-beta)*phiv.oldTime();
                const surfaceScalarField phi1 = beta*phiv.oldTime();
                
                sEqn = fvScalarMatrix
                (
                    EulerDdt.fvmDdt(s_[is])
                  + fvc::div(phi0*upInterp(phi0, s_[is].oldTime()))
                  + upwindCon.fvmDiv(phi1, s_[is])
                );
                sEqn.solve();
                fluxCorr -= phi0*upInterp(phi0, s_[is].oldTime()) + sEqn.flux();
            }
            else
            {
                const surfaceScalarField phi0
                    = (1-beta)*phiv.oldTime()*upInterp(phiv.oldTime(),s_[0].oldTime());
                const surfaceScalarField phi1
                    = beta*phiv.oldTime()*upInterp(phiv.oldTime(), s_[0]);
                
                sEqn = fvScalarMatrix
                (
                    EulerDdt.fvmDdt(s_[0], s_[is])
                  + fvc::div(phi0*upInterp(phi0, s_[is].oldTime()))
                  + upwindCon.fvmDiv(phi1, s_[is])
                );
                sEqn.solve();
                fluxCorr -= phi0*upInterp(phi0, s_[is].oldTime()) + sEqn.flux();
            }
            
            //fvc::fluxLimit(s_[is], fluxCorr, dt, FCTiter_[is]);
        }
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
