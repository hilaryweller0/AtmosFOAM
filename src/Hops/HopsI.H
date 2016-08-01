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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::tmp<Foam::surfaceScalarField> Foam::Hops::ddirToFlux
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("Hops::ddirToFlux") << " called for field " << V.name()
            << " which has a different mesh to Hops"
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
            ddirToFluxOffDiag(V)
        )
    );
    
    // add central ddir component
    tU.ref() += V*Hdiag();
    
    return tU;
}

inline Foam::tmp<Foam::surfaceScalarField> Foam::Hops::ddirToFluxOffDiag
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("Hops::ddirToFluxOffDiag") << " called for field " << V.name()
            << " which has a different mesh to Hops"
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
    
    surfaceVectorField delHat = delta()/mag(delta());
    
    volSymmTensorField T(fvc::surfaceSum(sqr(delHat)*mesh_.magSf()));
    if (mesh_.nSolutionD() != 3)
    {
        T += sqr(vector(mesh_.solutionD()) - vector(1,1,1))*magSqr(mesh_.C());
    }
    
    surfaceVectorField uf
    (
        fvc::interpolate
        (
            inv(T)
          & (fvc::surfaceSum(delHat*V*mesh_.magSf()/mag(delta()))), "H"
        )
    );
    
    // Take dot product in Sf direction and remove central delta component
    U = uf & (mesh_.Sf() - delta()*Hdiag());
    
    return tU;
}

inline Foam::tmp<Foam::surfaceScalarField> Foam::Hops::ddirToFluxCorr
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("Hops::ddirToFluxCorr") << " called for field " << V.name()
            << " which has a different mesh to Hops"
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
    
    return tU;
}

inline Foam::tmp<Foam::volVectorField> Foam::Hops::reconstructd
(
    const surfaceScalarField& V
) const
{
    if (V.mesh() != mesh_)
    {
        FatalErrorIn("Hops::reconstructd") << " called for field "
            << V.name()
            << " which has a different mesh to Hops"
            << exit(FatalError);
    }

    volSymmTensorField T(fvc::surfaceSum(sqr(delta())));
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
                "reconstructd(" + V.name() + ')',
                V.instance(),
                V.db()
            ),
            inv(T)
          & (fvc::surfaceSum(delta()*V))
        )
    );
    
    return tu;    
}

inline Foam::tmp<Foam::volVectorField> Foam::Hops::reconstruct
(
    const surfaceScalarField& phi
) const
{
    if (phi.mesh() != mesh_)
    {
        FatalErrorIn("Hops::reconstruct") << " called for field "
            << phi.name()
            << " which has a different mesh to Hops"
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
