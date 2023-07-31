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

#include "fvFCTadvectionDiffusion.H"
#include "fvcFluxLimit.H"
#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "upwind.H"
#include "linear.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "correctedSnGrad.H"
#include "gaussDivScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

enum class FCTType{FCTon, FCToff, FCTfinal};
NamedEnum<FCTType, 3> FCT_Type;
template<> const char* Foam::NamedEnum<FCTType,3>::names[]
    = {"FCTon", "FCToff", "FCTfinal"};

void FCTadvectionDiffusion
(
    volScalarField& T,
    const volScalarField& rho,
    const surfaceScalarField& offCentre,
    const surfaceScalarField& flux,
    const surfaceScalarField& gamma,
    const volScalarField& S,
    const bool finalIter
)
{
    // Reference to the mesh and the time step
    const fvMesh& mesh = T.mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Input
    dictionary dict
    (
        dictionary
        (
            mesh.schemes().div("advectDiffuse("+flux.name()+','+T.name()+')')
        ).subDict("advectionDiffusion")
    );
    
    const label nIter(readLabel(dict.lookup("nIterations")));
    
    const FCTType FCT(FCT_Type.read(dict.lookup("FCTlimit")));
    const scalar FCTmin
    (
        FCT != FCTType::FCToff ? readScalar(dict.lookup("FCTmin")) : scalar(0)
    );
    const scalar FCTmax
    (
        FCT != FCTType::FCToff ? readScalar(dict.lookup("FCTmax")) : scalar(0)
    );
    
    ITstream corrSchemeIS(dict.lookup("correctionScheme"));
    OStringStream oss;
    oss << corrSchemeIS[0] << " " << corrSchemeIS[1] << " " << corrSchemeIS[2]
        << " " << corrSchemeIS[3].wordToken()+"_0";
    IStringStream corrSchemeIS0(oss.str());
    
    // Old and new fluxes
    surfaceScalarField aphi = offCentre*flux;
    surfaceScalarField bphi = (1-offCentre)*flux.oldTime();
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Schemes needed
    linear<scalar> linearInterp(mesh);
    EulerDdtScheme<scalar> backwardEuler(mesh);
    gaussConvectionScheme<scalar> upwindConvect
    (
        mesh,
        aphi,
        tmp<surfaceInterpolationScheme<scalar>>(new upwind<scalar>(mesh, aphi))
    );
    gaussConvectionScheme<scalar> upwindConvectOld
    (
        mesh,
        bphi,
        tmp<surfaceInterpolationScheme<scalar>>
        (
            new upwind<scalar>(mesh, bphi)
        )
    );
    tmp<snGradScheme<scalar>> snGrad
    (
        snGradScheme<scalar>::New(mesh, mesh.schemes().snGrad(T.name()))
    );
    gaussLaplacianScheme<scalar,scalar> laplacian
    (
        mesh, linearInterp, snGrad
    );
    
    // The HO divergence scheme
    gaussConvectionScheme<scalar> gaussScheme
    (
        mesh, aphi, corrSchemeIS
    );
    gaussConvectionScheme<scalar> gaussSchemeOld
    (
        mesh, bphi, corrSchemeIS0
    );
    
    // The correction scheme
    const tmp<surfaceInterpolationScheme<scalar>>
        HOscheme(gaussScheme.interpScheme());
    const tmp<surfaceInterpolationScheme<scalar>> 
        HOschemeOld(gaussSchemeOld.interpScheme());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // offCentre on cell centres
    volScalarField offCentreC = fvc::average(1-offCentre);

    // Time averaged source
    volScalarField Smid = (1-offCentreC)*S.oldTime() + offCentreC*S;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Full solution without FCT

    // Initialise the high and low order total fluxes (including rho T)
    surfaceScalarField fluxLow = 
        bphi*upwindConvectOld.interpolate(bphi, T.oldTime())
      - (1-offCentre)*gamma.oldTime()*mesh.magSf()*snGrad->snGrad(T.oldTime());
    
    surfaceScalarField fluxHigh = bphi*HOschemeOld->correction(T.oldTime());
    
    // Create and the matrix equation without the HOC and without FCT
    fvMatrix<scalar> TEqn
    (
        backwardEuler.fvmDdt(rho, T)
      + upwindConvect.fvmDiv(aphi, T)
      - laplacian.fvmLaplacian(offCentre*gamma, T)
      + fvc::div(fluxLow)
      == Smid
    );

    // Iterations of full solution (without FCT)
    for(label iter = 0; iter < nIter; iter++)
    {
        solve(TEqn == -fvc::div(fluxHigh + aphi*HOscheme->correction(T)));
    }

    // Apply FCT if required
    if (FCT == FCTType::FCTon || (FCT == FCTType::FCTfinal && finalIter))
    {
        // Store full flux
        surfaceScalarField flux = fluxLow + fluxHigh
             +  aphi*HOscheme->correction(T) + TEqn.flux();

        // Calculate 1st-order (bounded) solution
        TEqn.solve();
        
        fluxLow += TEqn.flux();

        // Remove low order fluxes
        flux -= fluxLow;

        // Limit the solution
        Foam::fvc::fluxLimit(flux, rho*T, dt);

        // Calculate the limited solution
        T = 
        (
            rho.oldTime()*T.oldTime() + dt*
            (
               -fvc::div(flux + fluxLow)
              + Smid
            )
        )/rho;
        T.correctBoundaryConditions();

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
