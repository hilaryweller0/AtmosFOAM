/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "CLUST.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::CLUST<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "CLUST::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

    const surfaceVectorField& Uf = this->Uf_;

    //const GeometricField<Type,fvsPatchField,surfaceMesh> gradf=fvc::snGrad(vf);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > Grad = fvc::grad(vf);
    
    forAll(sfCorr, facei)
    {
        const scalar un = (Uf[facei] & mesh.Sf()[facei])/mesh.magSf()[facei];
        const scalar ub = mag(un)/max(SMALL, mag(Uf[facei]));
        label celli = un > 0 ? mesh.owner()[facei] : mesh.neighbour()[facei];
        label cellj = celli == mesh.owner()[facei] ? mesh.neighbour()[facei]
                                                    : mesh.owner()[facei];
        sfCorr[facei] = 0.5*b_*ub*
        (
            ((mesh.C()[cellj] - mesh.C()[celli]) & Grad[celli])
          + vf[celli] - vf[cellj]
        );
    }

    return tsfCorr;
}


namespace Foam
{
    //makeSurfaceInterpolationScheme(CLUST);
    makeSurfaceInterpolationTypeScheme(CLUST, scalar);
    makeSurfaceInterpolationTypeScheme(CLUST, vector);
}

// ************************************************************************* //
