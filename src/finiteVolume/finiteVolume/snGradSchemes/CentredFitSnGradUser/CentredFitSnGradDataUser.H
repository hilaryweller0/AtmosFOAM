/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::CentredFitSnGradDataUser

Description
    Data for centred fit snGrad schemes

SourceFiles
    CentredFitSnGradDataUser.C

\*---------------------------------------------------------------------------*/

#ifndef CentredFitSnGradDataUser_H
#define CentredFitSnGradDataUser_H

#include "FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedCentredCellToFaceStencil;

/*---------------------------------------------------------------------------*\
                    Class CentredFitSnGradDataUser Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class CentredFitSnGradDataUser
:
    public FitData
    <
        CentredFitSnGradDataUser<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
{
    // Private data

        //- For each cell in the mesh store the values which multiply the
        //  values of the stencil to obtain the gradient for each direction
        List<scalarList> coeffs_;


public:

    TypeName("CentredFitSnGradDataUser");


    // Constructors

        //- Construct from components
        CentredFitSnGradDataUser
        (
            const fvMesh& mesh,
            const extendedCentredCellToFaceStencil& stencil,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~CentredFitSnGradDataUser()
    {}


    // Member functions

        //- Return reference to fit coefficients
        const List<scalarList>& coeffs() const
        {
            return coeffs_;
        }

        //- Calculate the fit for the specified face and set the coefficients
        void calcFit
        (
            scalarList& coeffsi, // coefficients to be set
            const List<point>&,  // Stencil points
            const scalar wLin,   // Weight for linear approximation (weights
                                 // nearest neighbours)
            const scalar deltaCoeff, // uncorrected delta coefficient
            const label faci     // Current face index
        );

        void calcFit();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CentredFitSnGradDataUser.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
