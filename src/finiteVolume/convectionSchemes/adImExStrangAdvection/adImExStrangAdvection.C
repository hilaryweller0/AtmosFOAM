/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "adImExStrangAdvection.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "EulerDdtScheme.H"
#include "fvc.H"
//#include "fvcFluxLimit.H"
#include "localMax.H"
#include "CourantNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * * * //

template<class Type>
adImExStrangAdvection<Type>::adImExStrangAdvection
(
    const fvMesh& mesh,
    const surfaceScalarField& advFlux, // the advecting flux
    Istream& is
)
:
    convectionScheme<Type>(mesh, advFlux),
    dict_(is),
    RK_(dict_.lookup("RK_ButcherCoeffs")),
    gammaScale_(readScalar(dict_.lookup("gammaScale"))),
    gamma1cMax_(readScalar(dict_.lookup("gamma1cMax"))),
    gammaAdd_(gammaScale_ - gamma1cMax_),
    rhoName_(dict_.lookupOrDefault<word>("density", "none"))
    //alpha_(readScalar(dict_.lookup("alpha"))),
    //beta_(readScalar(dict_.lookup("beta"))),
    //gamma_(readScalar(dict_.lookup("gamma")))
    //fullSolver_(dict_.lookup("fullSolver")),
    //FCTlimit_(dict_.lookup("FCTlimit")),
    //FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
    //FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0)))
{}

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //


template<class Type>
tmp<SurfaceField<Type>>
adImExStrangAdvection<Type>::interpolate
(
    const surfaceScalarField& advFlux,
    const VolField<Type>& vf
) const
{
    // Un flux-corrected interpolation
    tmp<SurfaceField<Type>> tinterp
    (
        new SurfaceField<Type>
        (
            corrScheme(advFlux)->interpolate(vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<SurfaceField<Type>>
adImExStrangAdvection<Type>::flux
(
    const surfaceScalarField& advFlux,
    const VolField<Type>& vf
) const
{
    return advFlux*interpolate(advFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
adImExStrangAdvection<Type>::fvmDiv
(
    const surfaceScalarField& advFlux,
    const VolField<Type>& vf
) const
{
    // Matrix to be returned will only have a source term
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            advFlux.dimensions()*vf.dimensions()
        )
    );
    tfvm.ref() += fvcDiv(advFlux, vf);
    
    return tfvm;
}


template<class Type>
tmp<VolField<Type>>
adImExStrangAdvection<Type>::fvcDiv
(
    const surfaceScalarField& advFlux,
    const VolField<Type>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Schemes needed
    //localMax<scalar> maxInterp(mesh);
    EulerDdtScheme<Type> backwardEuler(mesh);
    
    // Create temporary field to advect
    VolField<Type> T
    (
        IOobject(vf.name(), mesh),
        vf,
        vf.boundaryField().types()
    );
    T.oldTime();
    
    // Calculate the Courant number on faces
    localMax<scalar> maxInterp(mesh);
    volScalarField Co0 = CourantNo(advFlux.oldTime(), dt);
    volScalarField Co1 = CourantNo(advFlux, dt);
    surfaceScalarField Co = max
    (
        maxInterp.interpolate(Co0), maxInterp.interpolate(Co1)
    );
    // Calculate alpha, beta and gamma as a function of the max Courant number
    surfaceScalarField alpha = 1 - 1/max(scalar(2), Co);
    surfaceScalarField beta  = 1 - 1/max(scalar(1), Co);
    //surfaceScalarField gamma = 6.5/(max(Co,2.5)+scalar(4));
    surfaceScalarField gamma = gammaScale_/(gammaAdd_ + max(Co, gamma1cMax_));
    Info << "Co goes from " << min(Co).value() << " to " << max(Co).value() << nl
         << "alpha goes from " << min(alpha).value() << " to " << max(alpha).value() << nl
         << "beta goes from " << min(beta).value() << " to " << max(beta).value() << nl
         << "gamma goes from " << min(gamma).value() << " to " << max(gamma).value() 
         << endl;
    
    // Total fluxes for RK stages
    PtrList<SurfaceField<Type>> F(RK_.n()+1);
    
    // The zeroth flux is for the explicit low order update
    F.set
    (
        0,
        (1-alpha)*beta*advFlux.oldTime()*upInterp(advFlux.oldTime(),T)
    );

    // Accumulate the total flux, starting from the old time step terms
    surfaceScalarField totalFlux = F[0];
    
    if (rhoName_ == "none")
    {
        // Advance T by the old time step low order flux
        T -= dt*fvc::div(F[0]);

        // RK advection
        for(int iRK = 0; iRK < RK_.n(); iRK++)
        {
            // Sub-stage size
            scalar c = 0;
            for(int j = 0; j <= iRK; j++) c += RK_[iRK][j];
            
            // Advecting flux at the sub time
            const surfaceScalarField af = (1-c)*advFlux.oldTime() + c*advFlux;
            
            // New sub stage flux
            F.set(iRK+1, af*((1-beta)*upInterp(af, T) + gamma*hCorr(af, T)));
            
            // Update the flux sum
            totalFlux = F[0];
            for(int j = 0; j <= iRK; j++)
            {
                totalFlux += RK_[iRK][j]*F[j+1];
            }
        
            // Update T
            T = T.oldTime() - dt*fvc::div(totalFlux);
        }
        
        // Final implicit stage
        gaussConvectionScheme<Type> upwindCon(mesh, advFlux, upwindScheme(advFlux));
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(T)
          + fvc::div(totalFlux)
          + upwindCon.fvmDiv(alpha*beta*advFlux, T)
        );
        TEqn.solve();
        totalFlux += TEqn.flux();
    }
    else
    {
        const volScalarField& rho = mesh.lookupObjectRef<volScalarField>(rhoName_);
        VolField<Type> rhoT = rho.oldTime()*T.oldTime();
        
        // Advance rhoT by the old time step low order flux and keep T consistent
        rhoT -= dt*fvc::div(F[0]);
        T = rhoT/rho.oldTime();

        // RK advection
        for(int iRK = 0; iRK < RK_.n(); iRK++)
        {
            // Sub-stage size
            scalar c = 0;
            for(int j = 0; j <= iRK; j++) c += RK_[iRK][j];
            
            // Advecting flux at the sub time
            const surfaceScalarField af = (1-c)*advFlux.oldTime() + c*advFlux;
            
            // New sub stage flux
            F.set(iRK+1, af*((1-beta)*upInterp(af, T) + gamma*hCorr(af, T)));
            
            // Update the flux sum
            totalFlux = F[0];
            for(int j = 0; j <= iRK; j++)
            {
                totalFlux += RK_[iRK][j]*F[j+1];
            }
        
            // Update rhoT
            rhoT = rho.oldTime()*T.oldTime() - dt*fvc::div(totalFlux);
        }
        
        // Final implicit stage
        gaussConvectionScheme<Type> upwindCon(mesh, advFlux,upwindScheme(advFlux));
        fvMatrix<Type> TEqn
        (
            backwardEuler.fvmDdt(rho,T)
          + fvc::div(totalFlux)
          + upwindCon.fvmDiv(alpha*beta*advFlux, T)
        );
        TEqn.solve();
        totalFlux += TEqn.flux();
    }
    
    // Calculate the implied divergence and return
    return tmp<VolField<Type>>
    (
        new VolField<Type>
        (
            "div("+advFlux.name()+','+vf.name()+')',
            fvc::div(totalFlux)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
