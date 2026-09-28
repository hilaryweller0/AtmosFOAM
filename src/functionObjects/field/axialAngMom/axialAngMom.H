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
    Foam::functionObjects::axialAngMom

Description
    Evaluates and stores Shallow Water Eqn Energy

    Example:
    \verbatim
        #includeFunc axialAngMom()
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    axialAngMom.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_axialAngMom_H
#define functionObjects_axialAngMom_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFields.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class axialAngMom Declaration
\*---------------------------------------------------------------------------*/

class axialAngMom
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects

{

private:

    // Private Data

        //- The vector pointing east
        const volVectorField east_;

        //- The initial, total axial angular momentum
        dimensionedScalar axialAngMom0_;

        //- The axial component of the angular momentum
        volScalarField axialAngMom_;
        
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate axialAngMom
        void  calcAxialAngMom();


public:

    //- Runtime type information
    TypeName("axialAngMom");


    // Constructors

        //- Construct from Time and dictionary
        axialAngMom
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        axialAngMom(const axialAngMom&) = delete;


    //- Destructor
    virtual ~axialAngMom();


    // Member Functions

        //- Read the axialAngMom data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute at the start of the run
        virtual bool executeAtStart() const { return true; }

        //- Calculate the axialAngMom
        virtual bool execute();

        //- Write the axialAngMom
        virtual bool write();

    // Access
    
        //- The energy
        const volScalarField& field() const { return axialAngMom_;}
        volScalarField& field() { return axialAngMom_;}

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const axialAngMom&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
