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

Namespace
    Foam

Class
    Foam::monitorFunction

Description
    Abstract base class for the monitor function for a Monge Ampere mesh
    from either a volScalarField defined on a different mesh or a specified
    function. Also provides a mapping from the monitor function onto a new
    set of points or cells

SourceFiles
    monitorFunction.C

\*---------------------------------------------------------------------------*/

#ifndef monitorFunction_H
#define monitorFunction_H

#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class monitorFunction Declaration
\*---------------------------------------------------------------------------*/

class monitorFunction
:
    public IOdictionary
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        monitorFunction(const monitorFunction&);

        //- Disallow default bitwise assignment
        void operator=(const monitorFunction&);


public:

    //- Runtime type information
    TypeName("monitorFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            monitorFunction,
            dictionary,
            (
                const IOdictionary& dict
            ),
            (dict)
        );


    // Constructors
    
        //- Construct from dictionary
        monitorFunction(const IOdictionary& dict);
    
    // Selectors
    
        //- Return a reference to the selected monitorFunction
        static autoPtr<monitorFunction> New
        (
            const IOdictionary& dict
        );

    //- Destructor
    virtual ~monitorFunction(){};
    
    // Member Functions
    
//        //- The scalarField of the monitor function at a Field of points
//        virtual tmp<scalarField> map
//        (
//            const pointField& pts,
//            const scalarField& oldMonitor
//        ) const = 0;
    
        //- The volScalarField of the monitor Function on the given mesh
        virtual tmp<volScalarField> map
        (
            const fvMesh& newMesh,
            const volScalarField& oldMonitor
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
