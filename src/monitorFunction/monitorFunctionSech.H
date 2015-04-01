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
    Foam::monitorFunctionSech

Description
    Defintes an analytic form of the monitor function using a sech function
    following Budd, Russell and Walsh, eqn 4.17

SourceFiles
    monitorFunctionSech.C

\*---------------------------------------------------------------------------*/

#ifndef monitorFunctionSech_H
#define monitorFunctionSech_H

#include "monitorFunction.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class monitorFunctionSech Declaration
\*---------------------------------------------------------------------------*/

class monitorFunctionSech
:
    public monitorFunction
{
    // coefficients of the Sech function
    scalar alpha1_, alpha2_, a_, sqra_;
    
    // centre of the radially symmetric feature
    point centre_;
    
public:

    //- Runtime type information
    TypeName("monitorFunctionSech");

    // Constructors

        //- Construct from components
        monitorFunctionSech
        (
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~monitorFunctionSech()
    {}


    // Member Functions

        // Access
        scalar alpha1() const { return alpha1_;}
        scalar alpha2() const { return alpha2_;}
        scalar a()  const { return a_;}
        scalar sqra()  const { return sqra_;}
        const point& centre() const { return centre_;}

//        //- The scalarField of the monitor function at a Field of points
//        tmp<scalarField> map
//        (
//            const pointField& pts,
//            const scalarField& oldMonitor
//        ) const;

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
