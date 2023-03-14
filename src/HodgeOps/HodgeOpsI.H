/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::tmp<Foam::surfaceScalarField> Foam::HodgeOps::ddirToFlux
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("HodgeOps::ddirToFlux") << " called for field " << V.name()
            << " which has a different mesh to HodgeOps"
            << exit(FatalError);
    }
    
    tmp<surfaceScalarField> tU
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddirToFlux(" + V.name() + ')',
                V.instance(),
                V.db()
            ),
            ddirToFluxOffDiag(V) + V*Hdiag()
        )
    );
    
    return tU;
}

inline Foam::tmp<Foam::surfaceScalarField> Foam::HodgeOps::ddirToFluxOffDiag
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("HodgeOps::ddirToFluxOffDiag") << " called for field " << V.name()
            << " which has a different mesh to HodgeOps"
            << exit(FatalError);
    }

    tmp<surfaceScalarField> tU
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddirToFlux(" + V.name() + ')',
                V.instance(),
                V.db()
            ),
            mesh_,
            V.dimensions()*dimLength
        )
    );
    surfaceScalarField& U = tU.ref();
    
    surfaceVectorField dhat = delta()/magd();
    //surfaceScalarField Sfdhat = mesh_.Sf() & dhat;
    
    //volTensorField T(fvc::surfaceSum(mesh_.Sf()*dhat));
    volSymmTensorField T(fvc::surfaceSum(sqr(dhat)*mesh_.magSf()));
    if (mesh_.nSolutionD() != 3)
    {
        T = T + sqr(vector(mesh_.solutionD()) - vector(1,1,1))
                *magSqr(mesh_.C());
    }
    
    surfaceVectorField uf
    (
        fvc::interpolate
        (
            inv(T)
//          & (fvc::surfaceSum(V*mesh_.Sf()/magd())), "H"
          & (fvc::surfaceSum(dhat*V*mesh_.magSf()/magd())), "H"
        )
    );
    
    // Take dot product in Sf direction and remove central delta component
    U = uf & (mesh_.Sf() - delta()*Hdiag());
    
    // Remove boundary components
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& pat = mesh_.boundary()[patchi];
        // only remove un-coupled boundary components
        if (!pat.coupled())
        {
            tU.ref().boundaryFieldRef()[patchi] == 0;
        }
    }

    return tU;
}

inline Foam::tmp<Foam::surfaceScalarField> Foam::HodgeOps::ddirToFluxCorr
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("HodgeOps::ddirToFluxCorr") << " called for field " << V.name()
            << " which has a different mesh to HodgeOps"
            << exit(FatalError);
    }

    tmp<surfaceScalarField> tU
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddirToFlux(" + V.name() + ')',
                V.instance(),
                V.db()
            ),
            ddirToFlux(V)
        )
    );
    
    // Remove uncorrected components
    tU.ref() -= V*mesh_.magSf()/mag(delta());
    
    // Remove boundary components
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& pat = mesh_.boundary()[patchi];
        // only remove un-coupled boundary components
        if (!pat.coupled())
        {
            tU.ref().boundaryFieldRef()[patchi] == 0;
        }
    }
    
    return tU;
}

inline Foam::tmp<Foam::volVectorField> Foam::HodgeOps::reconstructd
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("HodgeOps::reconstructd") << " called for field "
            << V.name()
            << " which has a different mesh to HodgeOps"
            << exit(FatalError);
    }

    surfaceVectorField dhat = delta()/magd();
    //surfaceScalarField Sfdhat = mesh_.Sf() & dhat;

    //volTensorField T(fvc::surfaceSum(mesh_.Sf()*dhat));
    volSymmTensorField T(fvc::surfaceSum(sqr(dhat)*mesh_.magSf()));
    if (mesh_.nSolutionD() != 3)
    {
        T = T
          + sqr(vector(mesh_.solutionD()) - vector(1,1,1))
               *magSqr(mesh_.C());
    }

    tmp<volVectorField> tu
    (
        new volVectorField
        (
            IOobject
            (
                "reconstructd(" + V.name() + ')',
                V.instance(),
                V.db()
            ),
            inv(T) & (fvc::surfaceSum(dhat*V*mesh_.magSf()/magd()))
        )
    );
    
    return tu;    
}

inline Foam::tmp<Foam::volVectorField> Foam::HodgeOps::reconstruct
(
    const surfaceScalarField& phi
) const
{
    if (phi.mesh() != mesh_)
    {
        FatalErrorIn("HodgeOps::reconstruct") << " called for field "
            << phi.name()
            << " which has a different mesh to HodgeOps"
            << exit(FatalError);
    }

    volSymmTensorField T(fvc::surfaceSum(sqr(mesh_.Sf())/mesh_.magSf()));
    if (mesh_.nSolutionD() != 3)
    {
        T += sqr(vector(mesh_.solutionD()) - vector(1,1,1))*magSqr(mesh_.C());
    }

    tmp<volVectorField> tu
    (
        new volVectorField
        (
            IOobject
            (
                "reconstruct(" + phi.name() + ')',
                phi.instance(),
                phi.db()
            ),
            inv(T)
          & (fvc::surfaceSum(mesh_.Sf()/mesh_.magSf()*phi)),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    
    tu.ref().correctBoundaryConditions();
    
    return tu;    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
