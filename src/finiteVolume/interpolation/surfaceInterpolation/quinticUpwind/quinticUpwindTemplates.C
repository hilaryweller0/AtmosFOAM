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

#include "quinticUpwind.H"
#include "orthogonalSnGrad.H"
//#include "localMax.H"
//#include "fvcLocalMinMax.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::quinticUpwind<Type>::quinticUpwind
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux
)
:
    upwind<Type>(mesh, faceFlux)
{}


template<class Type>
Foam::quinticUpwind<Type>::quinticUpwind
(
    const fvMesh& mesh,
    Istream& is
)
:
    upwind<Type>(mesh, is)
{}


template<class Type>
Foam::quinticUpwind<Type>::quinticUpwind
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    upwind<Type>(mesh, faceFlux, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::quinticUpwind<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tsfCorr
    (
        SurfaceField<Type>::New
        (
            "quinticUpwind::correction(" + vf.name() + ')',
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    SurfaceField<Type>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

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
    tmp<fv::gaussGrad<vector>> gradVector_
    (
        new fv::gaussGrad<vector>(mesh)
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        surfaceScalarField gradf = oSnGrad().snGrad
        (
            vf.component(cmpt), mesh.deltaCoeffs()
        );
        volVectorField gradV(gradScheme_().grad(vf.component(cmpt)));
        volTensorField gradGrad(gradVector_().grad(gradV));

        forAll(faceFlux, facei)
        {
            // Upwind and downwind cells
            const label u = (faceFlux[facei] >= 0) ? own[facei] : nei[facei];
            const label d = (faceFlux[facei] >  0) ? nei[facei] : own[facei];

            const vector dx = 2*(Cf[facei] - C[u]);
            const scalar magDx = dx & dhat[facei];

            setComponent(sfCorr[facei], cmpt) = 7/30.*gradf[facei]*magDx
                                              + 13/30.*(dx & gradV[u])
                                              - 1/6.*(dx & gradV[d])
                                              + 2/15.*(dx & (dx & gradGrad[u]));
        }
        
        // Boundary correction
        typename SurfaceField<Type>::
        Boundary& bSfCorr = sfCorr.boundaryFieldRef();
        
        forAll(bSfCorr, patchi)
        {
            fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];
            
            if (pSfCorr.coupled())
            {
                const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
                const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];
                
                const scalarField pgradf(gradf.boundaryField()[patchi]);
                
                const vectorField pGradVNei
                (
                    gradV.boundaryField()[patchi].patchNeighbourField()
                );
                const tensorField pGradGradNei
                (
                    gradGrad.boundaryField()[patchi].patchNeighbourField()
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
                        const vector dx = pd[facei];
                        const scalar magDx = mag(dx);
                        setComponent(pSfCorr[facei], cmpt) =
                            7/30.*pgradf[facei]*magDx
                          + 13/30.*(dx & gradV[own])
                          - 1/6.*(dx & pGradVNei[facei])
                          + 2/15.*(dx & (dx & gradGrad[own]));
                    }
                    else
                    {
                        const vector dx = -pd[facei];
                        const scalar magDx = -mag(dx);
                        setComponent(pSfCorr[facei], cmpt) =
                            7/30.*pgradf[facei]*magDx
                            + 13/30.*(dx & pGradVNei[facei])
                            - 1/6.*(dx & gradV[own])
                            + 2/15.*(dx & (dx & pGradGradNei[facei]));
                    }
                }
            }
        }
    }
    
    return tsfCorr;
}

// ************************************************************************* //
