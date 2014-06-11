/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::TVDlimited<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // references to some necessary data
    const fvMesh& mesh = vf.mesh();
    const vectorField& C = mesh.C();
    
    // Declare weights to return
    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                vf.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            scalar(0)
        )
    );
    surfaceScalarField& w = tw();

    // The unlimited interpolated field
    GeometricField<Type, fvsPatchField, surfaceMesh> sf
        = tInterp_().interpolate(vf);

    // Scalar version of the vol field
    volScalarField vsf = vf.component(0);
    // Scalar version of the interpolated field
    surfaceScalarField ssf = sf.component(0);
    
    // and its gradient
    fv::gaussGrad<scalar> grad(mesh);
    volVectorField gradvf(grad.grad(vsf));
//    typedef typename outerProduct<vector, scalar>::type GradType;
//    GeometricField<GradType, fvPatchField, volMesh> gradvf
//         = fv::gaussGrad<scalar>::gradf(sf, "name");
    
    //The TVDNVD limiter function
    NVDTVD tvd;
    
    // Calculate the weights for each (internal) sf
    forAll(ssf, face)
    {
        const label own = mesh.owner()[face];
        const label nei = mesh.neighbour()[face];
    
        // weight for central differencing
        const scalar cdWeight = mesh.surfaceInterpolation::weights()[face];
    
        // values needed by the limiter function
        const scalar faceFlux = faceFlux_[face];
        const scalar phiP = vsf[own];
        const scalar phiN = vsf[nei];
        const vector gradcP = gradvf[own];
        const vector gradcN = gradvf[nei];
        const vector& d = C[nei] - C[own];
        
        // the limiter function
        const scalar r = tvd.r(faceFlux, phiP, phiN, gradcP, gradcN, d);
        
        // Find the upwind value
        const scalar phiU = faceFlux > 0 ? phiP : phiN;
        
        // Face value from central differencing
        const scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;
        
        // Calculate the effective limiter for this interpolation
        scalar thisLimiter = stabilise(ssf[face] - phiU, SMALL)
                             /stabilise(phiCD - phiU, SMALL);

        // Limit the limiter to obey the TVD constraint
        thisLimiter = max(min(min(2*r, thisLimiter), 2), 0);
        //thisLimiter = max(min(thisLimiter, 2), 0);
        
//        Info << "limiter = " << thisLimiter << " r = " << r << nl;

        // Calculate the weight from the and the interpolation weight
        w[face] = thisLimiter*cdWeight
                + (1-thisLimiter)*pos(faceFlux_[face]);
        
//        Info << " cdWeight = " << cdWeight << " thisLimiter = "
//             << thisLimiter << " sf = " << sf[face] << " ssf = "
//             << ssf[face] << " w = " << w[face] << endl;
    }
    
    return tw;
}

// ************************************************************************* //
