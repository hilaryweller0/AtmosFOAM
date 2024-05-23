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
#include "EulerDdtScheme.H"
#include "upwind.H"
#include "linear.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "correctedSnGrad.H"
#include "gaussDivScheme.H"
#include "fvcAverage.H"

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


tmp<volScalarField> FCTadvectDiffuse
(
    const volScalarField& vf,
    const volScalarField& rho,
    const surfaceScalarField& fluxE,
    const surfaceScalarField& gamma,
    const volScalarField& S,
    int nIter,// = -1,
    const bool finalIter //= true
)
{
    // Reference to the mesh and the time step
    const fvMesh& mesh = vf.mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Create temporary field to advect
    volScalarField T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf,
        vf.boundaryField().types()
    );
    T.oldTime() == vf.oldTime();

    // Input
    dictionary dict
    (
        dictionary(mesh.schemes().div("advectDiffuse("+T.name()+')')).subDict("advectionDiffusion")
    );
    
    nIter = (nIter == -1) ?
        readLabel(dict.lookup("nIterations")):
        dict.lookupOrDefault<label>("nIterations", nIter);
    
    const FCTType FCT(FCT_Type.read(dict.lookup("FCTlimit")));
    const dimensionedScalar FCTmin
    (
        FCT != FCTType::FCToff ?
            dimensionedScalar(T.dimensions(), readScalar(dict.lookup("FCTmin"))) :
            dimensionedScalar(T.dimensions(), scalar(0))
    );
    const dimensionedScalar FCTmax
    (
        FCT != FCTType::FCToff ?
             dimensionedScalar(T.dimensions(), readScalar(dict.lookup("FCTmax"))) :
             dimensionedScalar(T.dimensions(), scalar(0))
    );
    
    ITstream corrSchemeIS(dict.lookup("correctionScheme"));
    
    // Schemes needed
    linear<scalar> linearInterp(mesh);
    gaussConvectionScheme<scalar> upwindConvect
    (
        mesh,
        fluxE,
        tmp<surfaceInterpolationScheme<scalar>>(new upwind<scalar>(mesh, fluxE))
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
    gaussConvectionScheme<scalar> gaussScheme(mesh, fluxE, corrSchemeIS);
    
    // The correction scheme
    const tmp<surfaceInterpolationScheme<scalar>>
        HOscheme(gaussScheme.interpScheme());

    // Store the low order flux
    surfaceScalarField fluxLow = fluxE*upwindConvect.interpolate(fluxE, T.oldTime())
                               - gamma*snGrad->snGrad(T.oldTime())*mesh.magSf();

    // Iterations of full solution (without FCT)
    surfaceScalarField fluxET = fluxE*HOscheme->correction(T);
    for(label iter = 0; iter < nIter; iter++)
    {
        if (iter > 0)
        {
            fluxET = fluxE*HOscheme->correction(T)
                   + fluxE*upwindConvect.interpolate(fluxE, T)
                   - fluxLow;
        }
        T = (rho.oldTime()*T.oldTime() - dt*
        (
            fvc::div(fluxLow + fluxET) - S
        ))/rho;
    }
    
    Info << "Before explicit FCT " << T.name() << " goes from " << min(T).value()
         << " to " << max(T).value() << endl;

    // Apply FCT if required
    if (FCT == FCTType::FCTon || (FCT == FCTType::FCTfinal && finalIter))
    {
        // Calculate 1st-order (bounded) solution
        T = (rho.oldTime()*T.oldTime() - dt*(fvc::div(fluxLow) - S))/rho;

        Info << "After explicit 1st order solution " << T.name() << " goes from "
             << min(T).value() << " to " << max(T).value() << endl;

        // Limit the solution
        if (FCTmin < FCTmax)
            {Foam::fvc::fluxLimit(fluxET, rho*T, rho*FCTmin, rho*FCTmax, dt);}
        else
            {Foam::fvc::fluxLimit(fluxET, rho*T, rho.oldTime()*T.oldTime(), dt);}

        // Testing (remove)
        T += -dt*fvc::div(fluxET)/rho;
        Info << "After explicit FCT " << T.name() << " goes from "
             << min(T).value() << " to " << max(T).value() << endl;
    }
    
    return tmp<volScalarField> 
    (
        new volScalarField
        (
            "ddt("+rho.name()+","+T.name()+")", 
            -fvc::div(fluxET + fluxLow) + S
        )
    );
}

tmp<volScalarField> FCTadvectDiffuse
(
    const volScalarField& vf,
    const volScalarField& rho,
    const surfaceScalarField& fluxI,
    const surfaceScalarField& fluxE,
    const surfaceScalarField& gamma,
    const volScalarField& S,
    int nIter,// = -1,
    const bool finalIter //= true
)
{
    // Reference to the mesh and the time step
    const fvMesh& mesh = vf.mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    // Create temporary field to advect
    volScalarField T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf,
        vf.boundaryField().types()
    );
    T.oldTime() == vf.oldTime();

    // The total flux
    const surfaceScalarField fluxSum = fluxI + fluxE;
    
    // Input
    dictionary dict
    (
        dictionary(mesh.schemes().div("advectDiffuse("+T.name()+')')).subDict("advectionDiffusion")
    );
    
    nIter = (nIter == -1) ?
        readLabel(dict.lookup("nIterations")):
        dict.lookupOrDefault<label>("nIterations", nIter);
    
    const FCTType FCT(FCT_Type.read(dict.lookup("FCTlimit")));
    const dimensionedScalar FCTmin
    (
        FCT != FCTType::FCToff ?
            dimensionedScalar(T.dimensions(), readScalar(dict.lookup("FCTmin"))) :
            dimensionedScalar(T.dimensions(), scalar(0))
    );
    const dimensionedScalar FCTmax
    (
        FCT != FCTType::FCToff ?
             dimensionedScalar(T.dimensions(), readScalar(dict.lookup("FCTmax"))) :
             dimensionedScalar(T.dimensions(), scalar(0))
    );
    
    ITstream corrSchemeIS(dict.lookup("correctionScheme"));
    
    // Schemes needed
    linear<scalar> linearInterp(mesh);
    EulerDdtScheme<scalar> backwardEuler(mesh);
    gaussConvectionScheme<scalar> upwindConvect
    (
        mesh,
        fluxSum,
        tmp<surfaceInterpolationScheme<scalar>>(new upwind<scalar>(mesh, fluxSum))
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
    gaussConvectionScheme<scalar> gaussScheme(mesh, fluxSum, corrSchemeIS);
    
    // The correction scheme
    const tmp<surfaceInterpolationScheme<scalar>>
        HOscheme(gaussScheme.interpScheme());

    // Store the low order flux
    surfaceScalarField fluxLow = fluxE
        *upwindConvect.interpolate(fluxE, T.oldTime());

    // Create and the matrix equation without the HOC and without FCT
    fvMatrix<scalar> TEqn
    (
        backwardEuler.fvmDdt(rho, T)
      + upwindConvect.fvmDiv(fluxI, T)
      + fvc::div(fluxLow)
      - laplacian.fvmLaplacian(gamma, T)
      - S
    );

    // Iterations of full solution (without FCT)
    surfaceScalarField fluxET = fluxSum*HOscheme->correction(T);
    for(label iter = 0; iter < nIter; iter++)
    {
        if (iter > 0)
        {
            fluxET = fluxSum*HOscheme->correction(T)
                   + fluxE*upwindConvect.interpolate(fluxE, T)
                   - fluxLow;
        }
        solve(TEqn == -fvc::div(fluxET));
    }
    
    Info << "Before FCT " << T.name() << "old goes from " << min(T.oldTime()).value()
         << " to " << max(T.oldTime()).value() << endl;
    Info << "Before FCT " << T.name() << " goes from " << min(T).value()
         << " to " << max(T).value() << endl;

    // Apply FCT if required
    if (FCT == FCTType::FCTon || (FCT == FCTType::FCTfinal && finalIter))
    {
        // Store full flux (without fluxLow)
        fluxET += TEqn.flux();
        
        // Add bits to T one at a time to see where it goes wrong
        T = (T.oldTime()*rho.oldTime() + dt*S)/rho.oldTime();
        Info << "After S " << T.name() << " goes from "
             << min(T).value() << " to " << max(T).value() << endl;
        T = (T*rho.oldTime() - dt*fvc::div(fluxLow))/rho;
        Info << "After fluxLow " << T.name() << " goes from "
             << min(T).value() << " to " << max(T).value() << endl;
        
        // Calculate 1st-order (bounded) solution and flux
        //T = T.oldTime();
        TEqn.solve();
        Info << "After 1st order solve " << T.name() << " goes from "
             << min(T).value() << " to " << max(T).value() << endl;
        
        fluxLow += TEqn.flux();
        // Remove implicit low order fluxes from fluxET
        fluxET -= TEqn.flux();

        // Limit the solution
        if (FCTmin < FCTmax)
            {Foam::fvc::fluxLimit(fluxET, rho*T, rho*FCTmin, rho*FCTmax, dt);}
        else
            {Foam::fvc::fluxLimit(fluxET, rho*T, rho.oldTime()*T.oldTime(), dt);}
    }
    else
    {
        fluxLow += TEqn.flux();
    }
    
    return tmp<volScalarField> 
    (
        new volScalarField
        (
            "ddt("+rho.name()+","+T.name()+")", 
            -fvc::div(fluxET + fluxLow) + S
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
