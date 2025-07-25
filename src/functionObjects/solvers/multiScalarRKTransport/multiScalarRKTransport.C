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
#include "velocityField.H"

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


// * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::functionObjects::multiScalarRKTransport::setFlux()
{
    // The mid time
    scalar vTime = time_.time().value() - 0.5*time_.deltaTValue();
    
    // Two options:
    //       set from velocityField function
    //       or look up
    // First consider the velocityField function
    if (velocityDictName_ != "")
    {
        IOdictionary dict
        (
            IOobject
            (
                velocityDictName_, time_.constant(), mesh_, IOobject::MUST_READ
            )
        );
        
        if (readBool(dict.lookup("timeVaryingWind")))
        {
            autoPtr<velocityField> v;
            v = velocityField::New(dict);
            v->applyTo(phiv_, vTime);
        }
    }
    
    // If the velocity flux already exits
    else if (mesh_.foundObject<surfaceScalarField>(phivName_))
    {
        const surfaceScalarField& phiv
             = mesh_.lookupObjectRef<surfaceScalarField>(phivName_);
        phiv_ = 0.5*(phiv.oldTime() + phiv);
    }
    // Else the velocity flux stays the same
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
    c_(RK_.n(), scalar(0)),
    correctionSchemeName_(dict.lookup("correctionScheme")),
    FCTiter_(dict.lookup("FCTiters")),
    FCTlimits_
    (
        dict.lookupOrDefault<scalarListList>
        (
            "FCTlimits",
            scalarListList(nFields_)
        )
    ),
    velocityDictName_
    (
        dict.lookupOrDefault<word>("velocityFieldDict", "")
    ),
    phiv_
    (
        mesh_.foundObject<surfaceScalarField>(phivName_) ?
        mesh_.lookupObjectRef<surfaceScalarField>(phivName_) :
        surfaceScalarField
        (
            IOobject
            (
                phivName_,
                runTime.startTime().name(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_
        )
    )
{
    if (nFields_ != fieldNames_.size() || nFields_ != FCTiter_.size())
    {
        FatalErrorIn("multiScalarRKTransport::multiScalarRKTransport")
            << "nFields = " << nFields_ << " but fields are " << fieldNames_
            << " and FCTiters are " << FCTiter_
            << exit(FatalError);
    }
    
    for(label iRK = 0; iRK < RK_.n(); iRK++)
    {
        for(label j = 0; j <= iRK; j++)
        c_[iRK] += RK_[iRK][j];
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
    // Update the velocity if needed
    setFlux();
    
    // Schemes needed
    localMax<scalar> maxInterp(mesh_);
    fv::EulerDdtScheme<scalar> EulerDdt(mesh_);
    fv::gaussConvectionScheme<scalar> upwindCon(mesh_,phiv_,upwindScheme(phiv_));
    
    // Reference to the time step
    const dimensionedScalar& dt = time_.deltaT();

    // Calculate the maximum Courant number on faces
    surfaceScalarField Co = maxInterp.interpolate(CourantNo(phiv_, dt));
    // Calculate alpha as a function of the max Courant number
    surfaceScalarField alpha = 1 - 1/max(scalar(2), Co);
    surfaceScalarField beta = 1 - 1/max(scalar(1), Co);
    surfaceScalarField gamma = 1 - 1/max(scalar(GREAT), Co);
    Info << "Co goes from " << min(Co).value() << " to " << max(Co).value() <<nl
        <<"alpha goes from "<<min(alpha).value()<<" to "<<max(alpha).value()<<nl
        <<"beta goes from "<<min(beta).value()<<" to "<<max(beta).value()<<nl
     <<"gamma goes from "<<min(gamma).value()<<" to "<<max(gamma).value()<<endl;

    // Store the high and low-order face values for each stage, for each field
     PtrList<PtrList<surfaceScalarField>> sL(RK_.n());
     PtrList<PtrList<surfaceScalarField>> sHC(RK_.n());
     PtrList<surfaceScalarField> rhoX(RK_.n());
         for(label iRK = 0; iRK < RK_.n(); iRK++)
    {
        sL.set(iRK, new PtrList<surfaceScalarField>(nFields_));
        sHC.set(iRK, new PtrList<surfaceScalarField>(nFields_));
    }

    // RK advection stages
    for(int iRK = 0; iRK < RK_.n(); iRK++)
    {
        // Update each tracer field, s_
        for(int is = 0; is < nFields_; is++)
        {
            // Is this field density weighted?
            const bool densityWeighted = is > 0 && massFluxName_ != "none";

            // Calculate the high and low-order face values
            sL[iRK].set(is, upInterp(phiv_, s_[is]));
            sHC[iRK].set(is, hCorr(phiv_, s_[is]));

            // Sum the explicit fluxes
            surfaceScalarField F = c_[iRK]*(1-alpha)*beta*sL[0][is];
            
            if (!densityWeighted)
            {
                for(label j = 0; j <= iRK; j++)
                {
                    F += RK_[iRK][j]*((1-beta)*sL[j][is] + gamma*sHC[j][is]);
                }
            }
            else
            {
                F *= rhoX[0];
                for(label j = 0; j <= iRK; j++)
                {
                    F += RK_[iRK][j]*rhoX[j]*
                        ((1-beta)*sL[j][is] + gamma*sHC[j][is]);
                }
            }
            
            // Assemble and solve the transport equation
            volScalarField divF("divF", fvc::div(phiv_*F));
            if (!densityWeighted)
            {
                s_[is] = s_[is].oldTime() - dt*divF;
                s_[is][0] += 1;
                s_[is][1] -= 1;
                fvScalarMatrix sEqn
                (
                    EulerDdt.fvmDdt(s_[is])
                  + divF
                  + upwindCon.fvmDiv(c_[iRK]*alpha*beta*phiv_, s_[is])
                );
                sEqn.solve();
            }
            else
            {
                s_[is] = (s_[0].oldTime()*s_[is].oldTime() - dt*divF)/s_[0];
                s_[is][0] += 1;
                s_[is][1] -= 1;
                fvScalarMatrix sEqn
                (
                    EulerDdt.fvmDdt(s_[0], s_[is])
                  + divF
                  + upwindCon.fvmDiv(c_[iRK]*alpha*beta*sL[iRK][0]*phiv_, s_[is])
                );
                sEqn.solve();
            }
            
            // Set the density face values for fluxes
            if (is == 0 && massFluxName_ != "none")
            {
                rhoX.set(0, sL[0][0] + RK_[iRK][0]*gamma*sHC[0][0]
                              /(c_[iRK]*beta*(1-alpha) + RK_[iRK][0]*(1-beta)));
                for(label j = 1; j <= iRK; j++)
                {
                    rhoX.set(j, sL[j][0] + gamma/(1-beta)*sHC[j][0]);
                }
            }
        }
    }

/*    // Accumulate the total flux for each scalar
    PtrList<surfaceScalarField> totalFlux(nFields_);
    for(int is = 0; is < nFields_; is++)
    {
        totalFlux.set(is, F[0][is]);
        for(label j = 0; j < RK_.n(); j++)
        {
            totalFlux[is] += RK_[RK_.n()-1][j]*F[j+1][is];
        }
    }
*/
    // Apply FCT if needed
/*    for(label is = 0; is < nFields_; is++)
    {
        if (FCTiter_[is] > 0)
        {
            surfaceScalarField fluxCorr = totalFlux[is];

            // Find the bounded solution and the bounded flux
            if (!densityWeighted)
            {
                const surfaceScalarField phi0 = (1-beta)*phiv_.oldTime();
                const surfaceScalarField phi1 = beta*phiv_.oldTime();

                volScalarField sInc
                     = fvc::div(phi0*upInterp(phi0, s_[is].oldTime()));
                s_[is] = s_[is].oldTime() - dt*sInc;
                     
           }
            else
            {
                volScalarField betaC = max
                (
                    1- s_[0].oldTime()/max
                    (
                        dt*fvc::surfaceIntegrateOut(totalFlux[0]),
                        dimensionedScalar(s_[0].dimensions(), SMALL)
                    ),
                    scalar(0)
                );
                beta = maxInterp.interpolate(betaC);
                Info << "New beta goes from " << min(beta).value() << " to "
                     << max(beta).value() << endl;

                surfaceScalarField phi0 = (1-beta)*totalFlux[0];
                surfaceScalarField phi1 = totalFlux[0] - phi0;

                volScalarField sInc
                    = fvc::div(phi0*upInterp(phi0, s_[is].oldTime()));
                s_[is] = (s_[0].oldTime()*s_[is].oldTime() - dt*sInc)/s_[0];
                
                fvScalarMatrix sEqn
                (
                    EulerDdt.fvmDdt(s_[0], s_[is])
                  + sInc
                  + upwindCon.fvmDiv(phi1, s_[is])
                );
                sEqn.solve();

                fluxCorr -= phi0*upInterp(phi0, s_[is].oldTime())
                          + sEqn.flux();
            }

            if (FCTlimits_[is].size() == 0)
            {
                fvc::fluxLimit(s_[is], fluxCorr, dt, FCTiter_[is]);
            }
            else if (FCTlimits_[is].size() == 1)
            {
                fvc::fluxLimit
                (
                    s_[is], fluxCorr, FCTlimits_[is][0], GREAT, dt, FCTiter_[is]
                );
            }
            else
            {
                fvc::fluxLimit
                (
                    s_[is], fluxCorr, FCTlimits_[is][0], FCTlimits_[is][1],
                    dt, FCTiter_[is]
                );
            }
        }
    }
*/
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
