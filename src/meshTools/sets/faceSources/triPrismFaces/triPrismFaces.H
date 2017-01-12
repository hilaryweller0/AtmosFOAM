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

Class
    Foam::triPrismFaces

Description
    A topoSetSource to select faces if owner and neighbour are both
    triangular prisms

SourceFiles
    triPrismFaces.C

\*---------------------------------------------------------------------------*/

#ifndef triPrismFaces_H
#define triPrismFaces_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class triPrismFaces Declaration
\*---------------------------------------------------------------------------*/

class triPrismFaces
:
    public topoSetSource
{

private:

        //- Add usage string
        static addToUsageTable usage_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("triPrismFaces");

    // Constructors

        //- Construct from components
        triPrismFaces
        (
            const polyMesh& mesh
        );

        //- Construct from dictionary
        triPrismFaces(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        triPrismFaces(const polyMesh& mesh, Istream&);


    //- Destructor
    virtual ~triPrismFaces();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACESETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
