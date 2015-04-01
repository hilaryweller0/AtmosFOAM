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
    Foam::monitorFunctionTanh

Description
    Defintes an analytic form of the monitor function using a tanh function
    following Ju and Gunzburger, MWR, 2011, eqn 4

SourceFiles
    monitorFunctionTanh.C

\*---------------------------------------------------------------------------*/

#ifndef monitorFunctionTanh_H
#define monitorFunctionTanh_H

#include "monitorFunction.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class monitorFunctionTanh Declaration
\*---------------------------------------------------------------------------*/

class monitorFunctionTanh
:
    public monitorFunction
{
    // coefficients of the tanh function
    scalar alpha_, beta_, gamma_;
    
    // centre of the high resoluiton
    point centre_;
    
    // specify spherical or Cartesian distance from the centre
    Switch spherical_;
    
public:

    //- Runtime type information
    TypeName("monitorFunctionTanh");

    // Constructors

        //- Construct from components
        monitorFunctionTanh
        (
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~monitorFunctionTanh()
    {}


    // Member Functions

        // Access
        scalar alpha() const { return alpha_;}
        scalar beta()  const { return beta_;}
        scalar gamma() const { return gamma_;}
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
