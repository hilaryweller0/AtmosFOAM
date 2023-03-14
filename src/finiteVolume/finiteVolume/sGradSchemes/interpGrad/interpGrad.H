/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::interpGrad

Description
    A surface grad scheme that simply interpolates the cell-centred
    gradient onto faces using interpolation and gradient schemes
    specified by fvSchemes.

SourceFiles
    interpGrad.C

\*---------------------------------------------------------------------------*/

#ifndef interpGrad_H
#define interpGrad_H

#include "sGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

template<class Type>
class interpGrad
:
    public fv::sGradScheme<Type>
{
public:
    TypeName("interpGrad");

    interpGrad(const fvMesh& mesh)
    :
        sGradScheme<Type>(mesh)
    {}

    interpGrad(const fvMesh& mesh, Istream&)
    :
        sGradScheme<Type>(mesh)
    {}

    virtual tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvsPatchField, surfaceMesh>
    > operator()
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    ) const;
    
private:
    interpGrad(const interpGrad&);
    void operator=(const interpGrad&);
};

}

}

#ifdef NoRepository
    #include "interpGrad.C"
#endif

#endif
