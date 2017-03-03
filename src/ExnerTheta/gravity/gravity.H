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
    Foam::gravity

Description
    A representation of gravitational acceleration on limited-area
    meshes in which g is a uniform vector.

SourceFiles
    gravity.C

\*---------------------------------------------------------------------------*/

#ifndef gravity_H
#define gravity_H

#include "fvCFD.H"

namespace Foam
{

class gravity
{
public:
    gravity(const fvMesh& mesh);
    virtual ~gravity();

    dimensionedVector operator()() const;

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    const surfaceScalarField& unitFaceNormal() const
    {
        return unitFaceNormal_;
    }

private:
    const fvMesh& mesh_;
    const IOdictionary dict_;
    const dimensionedVector g_;
    const surfaceScalarField unitFaceNormal_;

    gravity(const gravity&);
    void operator=(const gravity&);

};

}

#endif
