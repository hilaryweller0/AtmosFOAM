/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    Foam::functionObjects::uStar

Description
    Evaluates and outputs turbulence uStar for models. Values written to
    time directories as field 'uStar' or 'uStar.\<phase\>'.

    Example of function object specification:
    \verbatim
    uStar1
    {
        type        uStar;
        libs        ("libfieldFunctionObjects2.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: uStar           | yes        |
        phase    | phase name                 | no         |
    \endtable

    Note:
        Writing field 'uStar' is done by default, but it can be overridden by
        defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    uStar.C

\*---------------------------------------------------------------------------*/

#ifndef uStar_functionObject_H
#define uStar_functionObject_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class momentumTransportModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class uStar Declaration
\*---------------------------------------------------------------------------*/

class uStar
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private data

        //- Optional phase name
        word phaseName_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate uStar
        tmp<volScalarField> calcuStar
        (
            const momentumTransportModel& turbModel
        );


public:

    //- Runtime type information
    TypeName("uStar");


    // Constructors

        //- Construct from Time and dictionary
        uStar
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        uStar(const uStar&) = delete;


    //- Destructor
    virtual ~uStar();


    // Member Functions

        //- Read the uStar data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the uStar field
        virtual bool execute();

        //- Write the uStar field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uStar&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
