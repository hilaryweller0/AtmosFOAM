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
    Foam::functionObjects::massTransport

Description
    Evolves a scalar transport equation for density or a density weighted field.

    - To specify the field name set the \c field entry
    - To employ the same numerical schemes as another field set
      the \c schemesField entry,
    - The \c diffusivity entry can be set to \c none, \c constant, \c viscosity
    - A constant diffusivity is specified with the \c D entry,
    - If a momentum transport model is available and the \c viscosity
      diffusivity option specified an effective diffusivity may be constructed
      from the laminar and turbulent viscosities using the diffusivity
      coefficients \c alphal and \c alphat:
      \verbatim
          D = alphal*nu + alphat*nut
      \endverbatim

    Example for density:
    \verbatim
        #includeFunc massTransport
        (
            rho,
            advectingFlux=phiv,
            scalarFlux=phi // save phi= rho*phiv as phi in the database
        )
    \endverbatim
    Example for a density weighted field:
    \verbatim
        #includeFunc massTransport
        (
            s,
            advectingFlux=phi,
            diffusivity=none,
            density=rho
        )
    \endverbatim

    with scheme specification:
        div(phiv,rho)     Gauss linear;

    and solver specification:
        tracer
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0;
        }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    massTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massTransport_H
#define functionObjects_massTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class massTransport Declaration
\*---------------------------------------------------------------------------*/

class massTransport
:
    public fvMeshFunctionObject
{
public:

        //- Enumeration defining the type of the diffusivity
        enum class diffusivityType
        {
            none,
            constant,
            viscosity
        };

        //- Diffusivity type names
        static const NamedEnum<diffusivityType, 3> diffusivityTypeNames_;


private:

    // Private Data

        //- Name of field to process
        word fieldName_;

        //- Name of advecting flux field (optional, defaults to phi)
        word phiName_;

        //- Name of calculated scalar flux field (optional, defaults to none)
        word fluxName_;

        //- Name of density field (optional, defaults to none)
        word rhoName_;

        //- The type of diffusivity
        diffusivityType diffusivity_;

        //- Constant diffusivity coefficient (optional)
        scalar D_;
        
        //- Number of outer corrector iterations (optional)
        int nCorr_;

        //- Laminar diffusivity coefficient (optional)
        scalar alphal_;

        //- Turbulent diffusivity coefficient (optional)
        scalar alphat_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- The scalar field
        volScalarField s_;


    // Private Member Functions

        //- Return the diffusivity field
        tmp<volScalarField> D() const;


public:

    //- Runtime type information
    TypeName("massTransport");


    // Constructors

        //- Construct from Time and dictionary
        massTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        massTransport(const massTransport&) = delete;


    //- Destructor
    virtual ~massTransport();


    // Member Functions

        //- Read the massTransport data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do not execute at the start of the run
        virtual bool executeAtStart() const
        {
            return false;
        }

        //- Calculate the massTransport
        virtual bool execute();

        //- Write the updated scalar field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
