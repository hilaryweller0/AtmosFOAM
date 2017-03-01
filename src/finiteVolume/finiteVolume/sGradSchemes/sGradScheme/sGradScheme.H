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
    Foam::fv::sGradScheme

Description
    Abstract base class for surface gradient schemes.

SourceFiles
    sGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef sGradScheme_H
#define sGradScheme_H

#include "tmp.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

namespace fv
{

template<class Type>
class sGradScheme
:
    public tmp<sGradScheme<Type>>::refCount
{
    const fvMesh& mesh_;

    sGradScheme(const sGradScheme&);
    void operator=(const sGradScheme&);

public:

    virtual const word& type() const = 0;

    declareRunTimeSelectionTable
    (
        tmp,
        sGradScheme,
        Istream,
        (const fvMesh& mesh, Istream& schemeData),
        (mesh, schemeData)
    );

    sGradScheme(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    static tmp<sGradScheme<Type>> New
    (
        const fvMesh& mesh,
        Istream& schemeData
    );

    virtual ~sGradScheme();

    const fvMesh& mesh() const
    {
        return mesh_;
    }
    
    virtual tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvsPatchField, surfaceMesh>
    > operator()
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const = 0;
};

}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFvsGradTypeScheme(SS, Type)                                         \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            sGradScheme<Type>::addIstreamConstructorToTable<SS<Type> >          \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }


#define makeFvsGradScheme(SS)                                                   \
                                                                               \
makeFvsGradTypeScheme(SS, scalar)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
