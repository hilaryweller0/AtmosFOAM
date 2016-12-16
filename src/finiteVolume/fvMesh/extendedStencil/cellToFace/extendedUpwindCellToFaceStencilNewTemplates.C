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

#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::extendedUpwindCellToFaceStencilNew::weightedSum
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<scalar> >& ownWeights,
    const List<List<scalar> >& neiWeights
) const
{
    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type> > ownFld;
    collectData(ownMap(), ownStencil(), fld, ownFld);
    List<List<Type> > neiFld;
    collectData(neiMap(), neiStencil(), fld, neiFld);

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>
            (
                fld.name(),
                fld.dimensions(),
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsfCorr.ref();

    // Internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        std::cerr << "flux across internalFace at " << mesh.time().timeName() << "\n";
        if (phi[faceI] >= 0)
        {
            if (mag(phi[faceI]) > SMALL) std::cerr << "phi >= 0" << "\n";
            // Flux out of owner. Use upwind (= owner side) stencil.
            const List<Type>& stField = ownFld[faceI];
            const List<scalar>& stWeight = ownWeights[faceI];

            forAll(stField, i)
            {
                sf[faceI] += stField[i]*stWeight[i];
            }
            if (mag(phi[faceI]) > SMALL) std::cerr << "corr " << sf[faceI] << "\n";
        }
        else
        {
            if (mag(phi[faceI]) > SMALL) std::cerr << "phi < 0" << "\n";
            const List<Type>& stField = neiFld[faceI];
            const List<scalar>& stWeight = neiWeights[faceI];

            forAll(stField, i)
            {
                sf[faceI] += stField[i]*stWeight[i];
            }
            if (mag(phi[faceI]) > SMALL) std::cerr << "corr " << sf[faceI] << "\n";
        }
    }

    // Boundaries. Either constrained or calculated so assign value
    // directly (instead of nicely using operator==)
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary bSfCorr = sf.boundaryField();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            std::cerr << "flux across " << pSfCorr.patch().name() << " at " << mesh.time().timeName() << "\n";
            label faceI = pSfCorr.patch().start();

            forAll(pSfCorr, i)
            {
                if (phi.boundaryField()[patchi][i] >= 0)
                {
                    std::cerr << "phi >= 0" << "\n";
                    // Flux out of owner. Use upwind (= owner side) stencil.
                    const List<Type>& stField = ownFld[faceI];
                    const List<scalar>& stWeight = ownWeights[faceI];

                    forAll(stField, j)
                    {
                        const scalar st = static_cast<const scalar>(stField[j]);
                        std::cerr << st << "*" << stWeight[j] << "\n";
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                    std::cerr << "corr " << pSfCorr[i] << "\n";
                }
                else
                {
                    std::cerr << "phi < 0" << "\n";
                    const List<Type>& stField = neiFld[faceI];
                    const List<scalar>& stWeight = neiWeights[faceI];

                    forAll(stField, j)
                    {
                        const scalar st = static_cast<const scalar>(stField[j]);
                        std::cerr << st << "*" << stWeight[j] << "\n";
                        pSfCorr[i] += stField[j]*stWeight[j];
                    }
                    std::cerr << "corr " << pSfCorr[i] << "\n";
                }
                faceI++;
            }
        }
    }

    return tsfCorr;
}


// ************************************************************************* //
