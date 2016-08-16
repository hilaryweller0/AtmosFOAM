/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::liquidPhase

Description
    Holds thermodynamic constants for a constant density liquid and a
    volScalarField of the volume fraction of the liquid

SourceFiles
    liquidPhase.C

\*---------------------------------------------------------------------------*/

#ifndef liquidPhase_H
#define liquidPhase_H

#include "ExnerTheta.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class liquidPhase Declaration
\*---------------------------------------------------------------------------*/

class liquidPhase
{
    // Private data

        //- Thermodynamic data
        const dimensionedScalar Cp_;
        const dimensionedScalar rho_;

        //- Concentration (volume fraction) field of the liquid
        volScalarField v_;
        
        //- Rate of change of liquid volume fraction
        volScalarField dvdt_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        liquidPhase(const liquidPhase&);

        //- Disallow default bitwise assignment
        void operator=(const liquidPhase&);


public:

    // Constructors

        //- Construct v from IOobject and mesh and construct other
        //- properties from dictionary
        liquidPhase
        (
            const IOobject&,
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    ~liquidPhase();


    // Member Functions

        // Access
            //- Additional thermodynamic data
            const dimensionedScalar& Cp() { return Cp_; }
            const dimensionedScalar& rho() { return rho_; }
            
            //- volume fraction
            volScalarField& v() { return v_; }
            //- Rate of change of liquid fraction
            volScalarField& dvdt() { return dvdt_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
