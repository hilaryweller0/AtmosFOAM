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

#include "LaxWendroff.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "fvc.H"
#include "CourantNoFunc.H"
#include "fvcFluxLimit.H"
#include "fvcLocalMinMax.H"
#include "EulerDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

template<class Type>
tmp<surfaceScalarField>
LaxWendroff<Type>::calcOffCentre
(
    const surfaceScalarField& volFlux
) const
{
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    
    tmp<surfaceScalarField> toff
    (
        new surfaceScalarField
        (
            IOobject("offCentre", mesh.time().timeName(), mesh),
            mesh,
            dimensionedScalar("", dimless, offCentre_)
        )
    );
    surfaceScalarField& off = toff.ref();
    
    if (offCentre_ < -SMALL)
    {
        // First set offCentre_ to Courant number (on the face)
        if (mesh.objectRegistry::template foundObject<volScalarField>("Co"))
        {
            off = maxInterp.interpolate
            (
                mesh.objectRegistry::template
                     lookupObject<volScalarField>("Co")
            );
        }
        else
        {
            off = maxInterp.interpolate(CourantNo(volFlux, dt));
        }
        // Next calculate off centering from Courant number
        off = max(scalar(0), 1 - 1/(off + 0.25));
    }
    
    return toff;
}


template<class Type>
void LaxWendroff<Type>::smoothFluxCorr
(
    surfaceScalarField& fluxCorr,
    const surfaceScalarField& offCentre
) const
{
    const fvMesh& mesh = this->mesh();

    surfaceVectorField V(linearInterpolate(fvc::reconstruct(fluxCorr)));
    surfaceScalarField imp = offCentre/(offCentre+SMALL);
    imp = maxInterp.interpolate(fvc::localMax(imp));
    fluxCorr = imp*(V & mesh.Sf()) + (1-imp)*fluxCorr;
}


template<class Type>
void LaxWendroff<Type>::calculateCorr
(
    surfaceScalarField& Tcorr,
    const surfaceScalarField& volFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& offCentre
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();
    const surfaceScalarField& rdelta = mesh.deltaCoeffs();

    // Schemes needed
    fv::orthogonalSnGrad<Type> snGrad(mesh);
    fv::gaussGrad<Type> grad(mesh);
    upwind<vector> upwindV(mesh, volFlux);
    //QUICKupwind<Type> HOcorr(mesh, volFlux);

    // Gradient of T in the cell centre to cell centre direction
    surfaceScalarField snGradT = snGrad.snGrad(vf);

    // Full face gradient of T: average between face and upwind location
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    > gradT = linearInterpolate(grad.grad(vf));
    gradT += (snGradT - (gradT & mesh.delta())*rdelta)
         * mesh.delta()*rdelta;
    gradT = 0.5*(gradT + upwindV.interpolate(grad.grad(vf)));
    
    // The correction in space
    //fluxCorr = 0.5*mag(faceFlux)*snGradT/rdelta;
    Tcorr = gradT & upwindDelta(volFlux);
    //Tcorr = HOcorr.correction(vf);
    
    // Apply a correction in time if needed
    if (timeCorrector_ == "advective")
    {
        // The full velocity field from the flux and correct
        surfaceVectorField Uf = linearInterpolate(fvc::reconstruct(volFlux));
        Uf += (volFlux - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

        surfaceScalarField uDotGradT(Uf & gradT);
        /*if (isRho())
        {
            uDotGradT += linearInterpolate(fvc::div(volFlux));
        }*/

        // add advetive form time correction
        Tcorr -= max(1-2*offCentre, scalar(0))*0.5*dt*uDotGradT;
    }
    else if(timeCorrector_ == "flux")
    {
        Tcorr -= max(1-2*offCentre, scalar(0))*0.5*dt
             *linearInterpolate(fvc::div(volFlux*linearInterpolate(vf)));
    }
}


template<class Type>
void LaxWendroff<Type>::advect
(
    volScalarField& T,
    const surfaceScalarField& faceFlux,
    const surfaceScalarField& volFlux
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();

    // Schemes needed
    EulerDdtScheme<Type> backwardEuler(mesh);
    
    // Calculate the local off centering for each face from the Courant number
    surfaceScalarField offCentre(calcOffCentre(volFlux));

    const Switch anyImplicit = max(offCentre).value() > SMALL;
    
    // Calculte the first-order part of the advection
    fvMatrix<Type> TEqn(T, faceFlux.dimensions()*T.dimensions());
    TEqn += GaussUpwind.fvcDiv((1-offCentre)*faceFlux, T.oldTime());
    if (isRho())
    {
        TEqn += fvMatrix<Type>(backwardEuler.fvmDdt(rho(), T));
    }
    else TEqn += fvMatrix<Type>(backwardEuler.fvmDdt(T));
    
    if (anyImplicit)
    {
        TEqn += fvMatrix<Type>(GaussUpwind.fvmDiv(offCentre*faceFlux, T));
        TEqn.solve();
    }
    else
    {
        T = TEqn.H()/TEqn.A();
    }
    
    // Calculate and apply the correction
    surfaceScalarField Tcorr
    (
        IOobject("Tcorr", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(T.dimensions(), scalar(0))
    );
    
    calculateCorr(Tcorr, volFlux, T, offCentre);

    // Make Tcorr into the fluxCorr
    Tcorr *= faceFlux;
    surfaceScalarField& fluxCorr(Tcorr);
    
    // Smooth fluxCorr where offCentre>0
    smoothFluxCorr(fluxCorr, offCentre);

    // Limit the fluxes with Zalesak FCT limiter
    if (FCTlimit_)
    {
        dimensionedScalar FCTmin("FCTmin", T.dimensions(), FCTmin_);
        dimensionedScalar FCTmax("FCTmax", T.dimensions(), FCTmax_);
        
        if (FCTmin_ < FCTmax_)
        {
            fvc::fluxLimit(fluxCorr, T, FCTmin, FCTmax, dt);
        }
        else if (!anyImplicit && !isRho())
        {
            fvc::fluxLimit(fluxCorr, T, T.oldTime(), dt);
        }
        else
        {
            fvc::fluxLimit(fluxCorr, T, dt);
        }
    }
    
    T -= dt*specific(fvc::div(fluxCorr));
}


template<class Type>
tmp<surfaceScalarField>
LaxWendroff<Type>::specific(const surfaceScalarField& flux) const
{
    if (!isRho())
    {
        return flux;
    }
    
    surfaceScalarField rhof(linearInterpolate(0.5*(rho() + rho().oldTime())));
    tmp<surfaceScalarField> tvolFlux(new surfaceScalarField(flux/rhof));
    return tvolFlux;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
LaxWendroff<Type>::specific
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (!isRho())
    {
        return vf;
    }

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>(vf/rho())
    );

    return tvf;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> 
LaxWendroff<Type>::specific
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    if (!isRho())
    {
        return tvf;
    }

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf2
    (
        new GeometricField<Type, fvPatchField, volMesh>(tvf()/rho())
    );
    tvf.clear();

    return tvf2;
}

template<class Type>
tmp<surfaceVectorField> LaxWendroff<Type>::upwindDelta
(
    const surfaceScalarField& flux
) const
{
    const fvMesh& mesh = flux.mesh();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<surfaceVectorField> tDelta
    (
        new surfaceVectorField
        (
            IOobject("upwindDelta", flux.instance(), mesh),
            mesh,
            dimLength
        )
    );
    surfaceVectorField& delta = tDelta.ref();
    
    forAll(flux, facei)
    {
        const label celli = (flux[facei] > 0) ? owner[facei] : neighbour[facei];
        delta[facei] = Cf[facei] - C[celli];
    }
    
    const GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary& bFlux
        = flux.boundaryField();
        
    forAll(bFlux, patchi)
    {
        const fvsPatchField<scalar>& pFlux = bFlux[patchi];
        fvsPatchField<vector>& pDelta = delta.boundaryFieldRef()[patchi];
        // d-vectors
        const vectorField pd(mesh.Cf().boundaryField()[patchi].patch().delta());
        
        if (pFlux.coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = Cf.boundaryField()[patchi];
            
            forAll(pOwner, facei)
            {
                label own = pOwner[facei];
                
                if (pFlux[facei] > 0)
                {
                    pDelta[facei] = pCf[facei] - C[own];
                }
                else
                {
                    pDelta[facei] = pCf[facei] - pd[facei] - C[own];
                }
            }
        }
    }
    
    return tDelta;
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
LaxWendroff<Type>::LaxWendroff
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    convectionScheme<Type>(mesh, faceFlux),
    dict_(is),
    offCentre_(readScalar(dict_.lookup("offCentre"))),
    timeCorrector_(dict_.lookupOrDefault<word>("timeCorrector", "flux")),
    FCTlimit_(dict_.lookup("FCTlimit")),
    FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
    FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0))),
    densityName_(dict_.lookupOrDefault<word>("density", "")),
    GaussUpwind
    (
        mesh,
        faceFlux,
        tmp<surfaceInterpolationScheme<Type>>
        (
            new upwind<Type>(mesh, faceFlux)
        )
    ),
    maxInterp(mesh)
{
    if (offCentre_ >= 0.5 - SMALL) timeCorrector_ = "none";

    else if
    (
        timeCorrector_ != "advective"
     && timeCorrector_ != "flux"
     && timeCorrector_ != "none"
    )
    {
        FatalErrorIn("LaxWendroff::LaxWendroff")
            << "timeCorrector must be one of advective, flux or none but "
            << timeCorrector_ << " was given" << exit(FatalError);
        
    }
    
    if (densityName_ == "" && faceFlux.dimensions()[0] > 0)
    {
        FatalErrorIn("LaxWendroff::LaxWendroff")
            << "no density name given but flux is a mass flux"
            << exit(FatalError);
    }
    
    if (densityName_ != "" && !isRho())
    {
        FatalErrorIn("LaxWendroff::LaxWendroff")
            << "density named " << densityName_ << " does not exist"
            << exit(FatalError);
    }
}

// * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
LaxWendroff<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Just use the low order interpolate
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            GaussUpwind.interpolate(faceFlux, vf)
        )
    );

    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
LaxWendroff<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
LaxWendroff<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Matrix to be returned will only have a source term
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    tfvm.ref() += fvcDiv(faceFlux, vf);
    
    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
LaxWendroff<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Reference to the mesh, time step and mesh spacing
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& dt = mesh.time().deltaT();

    // Initialise the divergence to return
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "convection(" + faceFlux.name() + ',' + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            vf.dimensions()*faceFlux.dimensions()/dimVolume
        )
    );

    // Create temporary field to advect
    GeometricField<Type, fvPatchField, volMesh> T
    (
        IOobject(vf.name(), mesh.time().timeName(), mesh),
        vf.oldTime(),
        vf.boundaryField().types()
    );

    if (isRho())
    {
        const surfaceScalarField volFlux = specific(faceFlux);
        advect(T, faceFlux, volFlux);
        tConvection.ref() = (rho().oldTime()*vf.oldTime() - rho()*T)/dt;
    }
    else
    {
        advect(T, faceFlux, faceFlux);
        tConvection.ref() = (vf.oldTime() - T)/dt;
    }

    return tConvection;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
