/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2025 OpenFOAM Foundation
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
    Foam::functionObjects::ground

Description
    Evaluates and surface fields at the ground

    Example:
    \verbatim
        #includeFunc ground(z0=0.1)
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    ground.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ground_H
#define functionObjects_ground_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class ground Declaration
\*---------------------------------------------------------------------------*/

class ground
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects

{

private:

    // Private Data

        //- The patch number of the ground
        const label gPatchi_;

        //- Surface roughness height
        volScalarField z0_;

        //- Friction velocity
        scalarField uStar_;

    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate ground
        bool  calcGround();


public:

    //- Runtime type information
    TypeName("ground");


    // Constructors

        //- Construct from Time and dictionary
        ground
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ground(const ground&) = delete;


    //- Destructor
    virtual ~ground();


    // Member Functions

        //- Read the ground data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute at the start of the run
        virtual bool executeAtStart() const { return true; }

        //- Calculate the ground fields
        virtual bool execute();

        //- Write the BL fields
        virtual bool write();

    // Access
    
        //- The surface roughness
        const scalarField& z0() { return z0_.boundaryField()[gPatchi_];}
        
        //- The friction velocity
        const scalarField& uStar() { return uStar_;}

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ground&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
