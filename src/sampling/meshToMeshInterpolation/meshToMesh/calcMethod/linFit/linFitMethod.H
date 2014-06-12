/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::linFitMethod

Description
    Mesh-to-mesh interpolation from cell centres to cell centres using a 
    linear fit polynomial

    Not volume conservative.

SourceFiles
    linFitMethod.C

\*---------------------------------------------------------------------------*/

#ifndef linFitMethod_H
#define linFitMethod_H

#include "mapNearestMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class linFitMethod Declaration
\*---------------------------------------------------------------------------*/

class linFitMethod
:
    public mapNearestMethod
{
protected:

    // Protected Member Functions

        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Disallow default bitwise copy construct
        linFitMethod(const linFitMethod&);

        //- Disallow default bitwise assignment
        void operator=(const linFitMethod&);


public:

    //- Run-time type information
    TypeName("linFit");

    //- Construct from source and target meshes
    linFitMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~linFitMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing and weights
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                labelListList& tgtToTgtAddr,
                scalarListList& tgtToTgtWght
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
