/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#ifndef HoperatorI_H
#define HoperatorI_H

#include "Hoperator.H"

namespace Foam
{
namespace Hops
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

inline tmp<surfaceScalarField> circToFluxOffDiag(const surfaceScalarField& vS)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();

    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject
            (
                "circToFlux" + vS.name(),
                vS.time().timeName(),
                mesh
            ),
            mesh,
            vS.dimensions()
        )
    );
    surfaceScalarField& uS = tuS();

    surfaceVectorField Uf
    (
        inv(fvc::interpolate(fvc::surfaceSum(dualMesh.magSf()*sqr(dualMesh.jdir())), "H"))
      & fvc::interpolate(fvc::surfaceSum(dualMesh.jdir()*vS), "H")
    );
    
    // Take dot product in Sf direction and correct and remove ddir component
    uS = dualMesh.dualMap(Uf) & (mesh.Sf() - mesh.ddir()*mesh.Hdiag()*mesh.magSf());
    
    // Uncorrect ddir component
    //uS = (dualMesh.dualMap(Uf) & mesh.Sf()) - mesh.Hdiag()*dualMesh.dualFluxMap(vS);
    
    return tuS;
}


inline tmp<surfaceScalarField> circToFlux(const surfaceScalarField& vS)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();

    return mesh.Hdiag()*dualMesh.dualFluxMap(vS) + circToFluxOffDiag(vS);
}


inline tmp<surfaceScalarField> ddirToFluxOffDiag
(
    const surfaceScalarField& vS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));

    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddirToFlux(" + vS.name() + ')',
                vS.instance(),
                vS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vS.mesh(),
            vS.dimensions()
        )
    );
    surfaceScalarField& uS = tuS();
    
    surfaceVectorField Uf
    (
//        inv(fvc::interpolate(fvc::surfaceSum(mesh.magSf()*sqr(mesh.ddir())), "H"))
//      & fvc::interpolate(fvc::surfaceSum(mesh.ddir()*vS), "H")
        fvc::interpolate
        (
            inv(fvc::surfaceSum(mesh.magSf()*sqr(mesh.ddir())))
          & (fvc::surfaceSum(mesh.ddir()*vS)), "H"
        )

    );
    
    // Take dot product in Sf direction and remove central ddir component
    uS = Uf & (mesh.Sf() - mesh.ddir()*mesh.Hdiag()*mesh.magSf());
    
    return tuS;
}

inline tmp<surfaceScalarField> ddirToFlux
(
    const surfaceScalarField& vS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));

    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddirToFlux(" + vS.name() + ')',
                vS.instance(),
                vS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ddirToFluxOffDiag(vS)
        )
    );
    surfaceScalarField& uS = tuS();
    
    // add central ddir component
    uS += vS*mesh.Hdiag();
    
    return tuS;
}


inline tmp<surfaceVectorField> faceReconstruct
(
    const surfaceScalarField& uS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    tmp<surfaceVectorField> tUf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "faceReconstruct(" + uS.name() + ')',
                uS.instance(),
                uS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate
            (
                inv(fvc::surfaceSum(mesh.magSf()*sqr(mesh.idir())))
              & fvc::surfaceSum(mesh.idir()*uS),
                "H"
            )
        )
    );
    
//    // correct part in direction idir
//    surfaceVectorField& Uf = tUf();
//    Uf += (uS/mesh.magSf() - (Uf & mesh.idir()))*mesh.idir();
    
    return tUf;
}


inline tmp<surfaceVectorField> faceReconstruct
(
    const surfaceScalarField& uS,
    const surfaceVectorField& dir
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    tmp<surfaceVectorField> tUf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "faceReconstruct(" + uS.name() + ')',
                uS.instance(),
                uS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate
            (
                inv(fvc::surfaceSum(mesh.magSf()*sqr(dir)))
              & fvc::surfaceSum(dir*uS),
                "H"
            )
        )
    );
    
//    // correct part in direction dir
//    surfaceVectorField& Uf = tUf();
//    Uf += (uS/mesh.magSf() - (Uf & dir))*dir;

    return tUf;
}

inline tmp<volVectorField> cellReconstruct
(
    const surfaceScalarField& uS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    tmp<volVectorField> tU
    (
        new volVectorField
        (
            IOobject
            (
                "cellReconstruct(" + uS.name() + ')',
                uS.instance(),
                uS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(fvc::surfaceSum(mesh.magSf()*sqr(mesh.idir())))
          & fvc::surfaceSum(mesh.idir()*uS)
        )
    );
    
    return tU;
}


inline tmp<volVectorField> cellReconstruct
(
    const surfaceScalarField& uS,
    const surfaceVectorField& dir
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    tmp<volVectorField> tU
    (
        new volVectorField
        (
            IOobject
            (
                "faceReconstruct(" + uS.name() + ')',
                uS.instance(),
                uS.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(fvc::surfaceSum(mesh.magSf()*sqr(dir)))
          & fvc::surfaceSum(dir*uS)
        )
    );
    
    return tU;
}

} // End namespace Hops
} // End namespace Foam
#endif

