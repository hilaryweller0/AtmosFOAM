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

#include "cubicUpwind.H"
#include "orthogonalSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::surfaceVectorField>
Foam::cubicUpwind<Foam::vector>::correction
(
    const volVectorField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceVectorField> tsfCorr
    (
        surfaceVectorField::New
        (
            "cubicUpwind::correction(" + vf.name() + ')',
            mesh,
            dimensioned<vector>(vf.name(), vf.dimensions(), Zero)
        )
    );

    surfaceVectorField& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField dhat = mesh.delta()/mag(mesh.delta());

    tmp<fv::orthogonalSnGrad<vector>> oSnGrad
    (
        new fv::orthogonalSnGrad<vector>(mesh)
    );
    tmp<fv::gaussGrad<vector>> gradScheme_
    (
        new fv::gaussGrad<vector>(mesh)
    );

    const volTensorField gradVf = gradScheme_().grad(vf);

    surfaceTensorField gradVff = linearInterpolate(gradVf);
    surfaceVectorField sdGrad = oSnGrad().snGrad(vf, mesh.deltaCoeffs());
    gradVff += (sdGrad - (gradVff & dhat))*dhat;

    forAll(faceFlux, facei)
    {
        const label celli =
            (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];
        
        sfCorr[facei] = 1./3.*(Cf[facei] - C[celli])
                         & (2*gradVf[celli] + gradVff[facei]);
    }


    typename surfaceVectorField::Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchVectorField& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = Cf.boundaryField()[patchi];
            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const tensorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );
            
            const tensorField pGradVff(gradVff.boundaryField()[patchi]);

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own])
                                   & (2*gradVf[own] + pGradVff[facei])/3.;
                }
                else
                {
                    pSfCorr[facei] = (pCf[facei] - pd[facei] - C[own])
                                   & (2*pGradVfNei[facei] + pGradVff[facei])/3.;
                }
            }
        }
    }

    return tsfCorr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makelimitedSurfaceInterpolationScheme(cubicUpwind)
}

// ************************************************************************* //
