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

#ifndef unitVectors_H
#define unitVectors_H

#include "GeometricScalarField.H"

#define TEMPLATE \
    template<class vector, template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class vector, template<class> class PatchField, class GeoMesh>
inline void unitVector
(
    GeometricField<vector, PatchField, GeoMesh>& guf,
    const GeometricField<vector, PatchField, GeoMesh>& gf
)
{
    guf = gf/max(mag(gf), SMALL*max(mag(gf)));
}

template<class vector, template<class> class PatchField, class GeoMesh>
inline tmp<GeometricField<vector, PatchField, GeoMesh> > unitVector
(
    const GeometricField<vector, PatchField, GeoMesh>& gf
)
{
    tmp<GeometricField<vector, PatchField, GeoMesh> > tGuf
    (
        new GeometricField<vector, PatchField, GeoMesh>
        (
            IOobject
            (
                "unitVector(" + gf.name() + ')',
                gf.instance(),
                gf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gf.mesh(),
            dimless
        )
    );
    
    unitVector(tGuf(), gf);
    
    return tGuf;
}

template<class vector, template<class> class PatchField, class GeoMesh>
inline tmp<GeometricField<vector, PatchField, GeoMesh> > unitVector
(
    const tmp<GeometricField<vector, PatchField, GeoMesh> >& tgf
)
{
    const GeometricField<vector, PatchField, GeoMesh>& gf = tgf();

    tmp<GeometricField<vector, PatchField, GeoMesh> > tGuf
    (
        new GeometricField<vector, PatchField, GeoMesh>
        (
            IOobject
            (
                "unitVector(" + gf.name() + ')',
                gf.instance(),
                gf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gf.mesh(),
            dimless
        )
    );
    
    unitVector(tGuf(), gf);
    
    tgf.clear();
    
    return tGuf;
}

inline tmp<vectorField> unitVector(const vectorField& gf)
{
    tmp<vectorField> tGuf(new vectorField(gf/mag(gf)));
    return tGuf;
}

inline tmp<vectorField> unitVector(const tmp<vectorField>& tgf)
{
    tmp<vectorField> tGuf
    (
        new vectorField
        (
            tgf()
            /max
            (
                mag(tgf()),
                SMALL*max(mag(tgf()))
            )
        )
    );
    tgf.clear();
    return tGuf;
}

inline tmp<vectorField> unitVector(const vectorIOField& gf)
{
    tmp<vectorField> tGuf(new vectorField(gf/max(mag(gf), SMALL*max(mag(gf)))));
    return tGuf;
}

inline tmp<vectorField> unitVector(const tmp<vectorIOField>& tgf)
{
    tmp<vectorField> tGuf
    (
        new vectorField(tgf()/max(mag(tgf()), SMALL*max(mag(tgf()))))
    );
    tgf.clear();
    return tGuf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#endif
// ************************************************************************* //
