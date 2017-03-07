/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::volInterpolationScheme

Description
    Abstract base class for vol interpolation schemes.

SourceFiles
    volInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef volInterpolationScheme_H
#define volInterpolationScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

template<class Type>
class volInterpolationScheme
:
    public tmp<volInterpolationScheme<Type>>::refCount
{
public:
    TypeName("volInterpolationScheme");

    declareRunTimeSelectionTable
    (
        tmp,
        volInterpolationScheme,
        Mesh,
        (
            const fvMesh& mesh,
            Istream& schemeData
        ),
        (mesh, schemeData)
    );

    volInterpolationScheme(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    static tmp<volInterpolationScheme<Type>> New
    (
        const fvMesh& mesh,
        Istream& schemeData
    );

    virtual ~volInterpolationScheme();

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>
    interpolate
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const = 0;

private:
    const fvMesh& mesh_;
    volInterpolationScheme(const volInterpolationScheme&);
    void operator=(const volInterpolationScheme&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeVolInterpolationTypeScheme(SS, Type)                           \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
volInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \

#define makeVolInterpolationScheme(SS)                                     \
                                                                               \
makeVolInterpolationTypeScheme(SS, scalar)                                 \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
