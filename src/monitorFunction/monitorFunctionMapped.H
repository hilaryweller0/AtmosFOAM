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
    Foam::monitorFunctionMapped

Description
    Maps the monitor function from one mesh to another using meshToMesh

SourceFiles
    monitorFunctionMapped.C

\*---------------------------------------------------------------------------*/

#ifndef monitorFunctionMapped_H
#define monitorFunctionMapped_H

#include "monitorFunction.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class monitorFunctionMapped Declaration
\*---------------------------------------------------------------------------*/

class monitorFunctionMapped
:
    public monitorFunction
{
    
public:

    //- Runtime type information
    TypeName("monitorFunctionMapped");

    // Constructors

        //- Construct from components
        monitorFunctionMapped
        (
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~monitorFunctionMapped()
    {}


    // Member Functions

        //- The volScalarField of the monitor Function on the given mesh
        tmp<volScalarField> map
        (
            const fvMesh& newMesh,
            const volScalarField& oldMonitor
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
