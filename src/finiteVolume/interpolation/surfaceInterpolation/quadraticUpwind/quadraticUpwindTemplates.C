/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "quadraticUpwind.H"
#include "orthogonalSnGrad.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::quadraticUpwind<Type>::quadraticUpwind
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux
)
:
    upwind<Type>(mesh, faceFlux),
    Co1_(-1),
    Co2_(-1),
    applyBlend_(false),
    rhoName_("")
{}


template<class Type>
Foam::quadraticUpwind<Type>::quadraticUpwind
(
    const fvMesh& mesh,
    Istream& is
)
:
    upwind<Type>(mesh, is),
    Co1_(readScalar(is)),
    Co2_(readScalar(is)),
    applyBlend_(Co2_ >= Co1_ && Co1_ >= 0),
    rhoName_("")
{}


template<class Type>
Foam::quadraticUpwind<Type>::quadraticUpwind
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    upwind<Type>(mesh, faceFlux, is),
    Co1_(readScalar(is)),
    Co2_(readScalar(is)),
    applyBlend_(Co2_ >= Co1_ && Co1_ >= 0),
    rhoName_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::quadraticUpwind<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "quadraticUpwind::correction(" + vf.name() + ')',
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    SurfaceField<Type>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField dhat = mesh.delta()/mag(mesh.delta());

    tmp<fv::orthogonalSnGrad<scalar>> oSnGrad
    (
        new fv::orthogonalSnGrad<scalar>(mesh)
    );
    tmp<fv::gaussGrad<scalar>> gradScheme_
    (
        new fv::gaussGrad<scalar>(mesh)
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        volVectorField gradVf(gradScheme_().grad(vf.component(cmpt)));
        surfaceVectorField gradVff = linearInterpolate(gradVf);
        surfaceScalarField sdGrad = oSnGrad().snGrad
        (
            vf.component(cmpt), mesh.deltaCoeffs()
        );
        gradVff += (sdGrad - (gradVff & dhat))*dhat;

        forAll(faceFlux, facei)
        {
            const label celli =
                (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];

            setComponent(sfCorr[facei], cmpt) =
                0.5*(Cf[facei] - C[celli]) & (gradVf[celli] + gradVff[facei]);
        }

        typename SurfaceField<Type>::
            Boundary& bSfCorr = sfCorr.boundaryFieldRef();

        forAll(bSfCorr, patchi)
        {
            fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

            if (pSfCorr.coupled())
            {
                const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
                const vectorField& pCf = Cf.boundaryField()[patchi];
                const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

                const vectorField pGradVfNei
                (
                    gradVf.boundaryField()[patchi].patchNeighbourField()
                );

                const vectorField pGradVff
                (
                    gradVff.boundaryField()[patchi]
                );

                // Build the d-vectors
                const vectorField pd
                (
                    Cf.boundaryField()[patchi].patch().delta()
                );

                forAll(pOwner, facei)
                {
                    label own = pOwner[facei];

                    if (pFaceFlux[facei] > 0)
                    {
                        setComponent(pSfCorr[facei], cmpt) =
                            0.5*(pCf[facei] - C[own])
                          & (gradVf[own] + pGradVff[facei]);
                    }
                    else
                    {
                        setComponent(pSfCorr[facei], cmpt) =
                            0.5*(pCf[facei] - pd[facei] - C[own])
                          & (pGradVfNei[facei] + pGradVff[facei]);
                    }
                }
            }
        }
    }
    //Scale the correction by the blending factor
    if (applyBlend_)
    {
        sfCorr *= blendingFactor();
    }

    return tsfCorr;
}

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::quadraticUpwind<Type>::blendingFactor() const
{
    const fvMesh& mesh = this->mesh();
    tmp<surfaceScalarField> tUflux = this->faceFlux_;
    
    if (this->faceFlux_.dimensions() == dimMassFlux)
    {
        /*// Currently assume that the density field
        // corresponding to the mass-flux is named "rho" or "rho_0"
        const word& fl = this->faceFlux_.name();
        const word rhoName = (fl.replace("_0", "") != fl)
             ? "rho_0" : "rho";
        
        //Info << "Flux is " << fl << " rho is " << rhoName << endl;
        */
        const volScalarField& rho =
            mesh.objectRegistry::template lookupObject<volScalarField>
            (rhoName_);

        tUflux = this->faceFlux_/linearInterpolate(rho);
    }
    else if (this->faceFlux_.dimensions() != dimVolumetricFlux)
    {
        FatalErrorInFunction
            << "dimensions of faceFlux are not correct"
            << exit(FatalError);
    }

    return surfaceScalarField::New
    (
        "BlendingFactor",
        scalar(1) - max
        (
            min
            (
                (
                    mesh.time().deltaT()*mesh.deltaCoeffs()
                   *mag(tUflux)/mesh.magSf()
                  - Co1_
                )/max(Co2_ - Co1_, SMALL),
                scalar(1)
            ),
            scalar(0)
        )
    );
}

// ************************************************************************* //
