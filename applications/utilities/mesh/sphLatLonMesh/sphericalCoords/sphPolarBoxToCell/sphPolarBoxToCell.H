/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::sphPolarBoxToCell

Description
    A topoSetSource to select cells based on cell centres inside lat/lon/radius box.

SourceFiles
    sphPolarBoxToCell.C

\*---------------------------------------------------------------------------*/

#ifndef sphPolarBoxToCell_H
#define sphPolarBoxToCell_H

#include "topoSetSource.H"
#include "polarPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sphPolarBoxToCell Declaration
\*---------------------------------------------------------------------------*/

class sphPolarBoxToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- bounding box.
        polarPoint BLcorner_;
        polarPoint TRcorner_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("sphPolarBoxToCell");

    // Constructors

        //- Construct from components
        sphPolarBoxToCell
        (
            const polyMesh& mesh,
            const polarPoint& BLcorner,
            const polarPoint& TRcorner
        );

        //- Construct from dictionary
        sphPolarBoxToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        sphPolarBoxToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~sphPolarBoxToCell();


    // Member Functions

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
