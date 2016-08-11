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
    Foam::perfectGasPhase

Description
    Holds thermodynamic constants for a perfect gas and a volScalarField of the
    density

SourceFiles
    perfectGasPhase.C

\*---------------------------------------------------------------------------*/

#ifndef perfectGasPhase_H
#define perfectGasPhase_H

#include "ExnerTheta.H"
#include "volFields.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class perfectGasPhase Declaration
\*---------------------------------------------------------------------------*/

class perfectGasPhase
:
    public constTransport<hConstThermo<perfectGas<specie> > >
{
    // Private data

        //- Additional thermodynamic data
        const dimensionedScalar R_;
        const dimensionedScalar T0_;
        const dimensionedScalar p0_;
        const dimensionedScalar Cp_;
        const dimensionedScalar Cv_;
        const scalar kappa_;

        //- Density field of the phase
        volScalarField rho_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        perfectGasPhase(const perfectGasPhase&);

        //- Disallow default bitwise assignment
        void operator=(const perfectGasPhase&);


public:

    // Constructors

        //- Construct rho from IOobject and mesh and construct constTransport
        //- from Istream
        perfectGasPhase
        (
            const IOobject&,
            const fvMesh&,
            Istream&,
            const dictionary&
        );


    //- Destructor
    ~perfectGasPhase();


    // Member Functions

        // Access
            //- Additional thermodynamic data
            const dimensionedScalar& R() { return R_; }
            const dimensionedScalar& T0() { return T0_; }
            const dimensionedScalar& p0() { return p0_; }
            const dimensionedScalar& Cp() { return Cp_; }
            const dimensionedScalar& Cv() { return Cv_; }
            scalar kappa() { return kappa_; }
            
            //- Density
            volScalarField& rho() { return rho_; }

            //- Partial pressure (from perfect gas eqn of state)
            tmp<volScalarField> partialPressure(const volScalarField& T);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
