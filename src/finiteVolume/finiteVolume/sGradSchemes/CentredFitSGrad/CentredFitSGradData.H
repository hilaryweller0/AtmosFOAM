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
    Foam::CentredFitSGradData

Description
    Data for the SGrad scheme

SourceFiles
    CentredFitSGradData.C

\*---------------------------------------------------------------------------*/

#ifndef CentredFitSGradData_H
#define CentredFitSGradData_H

#include "FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedCentredCellToFaceStencil;

/*---------------------------------------------------------------------------*\
                    Class CentredFitSGradData Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class CentredFitSGradData
:
    public FitData
    <
        CentredFitSGradData<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
{
    // Private data
    
        //- Store the face deltas for the mesh
        const surfaceVectorField delta_;

        //- For each cell in the mesh store the values which multiply the
        //  values of the stencil to obtain the gradient for each direction
        List<scalarList> coeffs_;

protected:

    //- Find the cell to cell direction (i) and j and k directions for face faci
        void findFaceDirs
        (
            vector& idir,        // value changed in return
            vector& jdir,        // value changed in return
            vector& kdir,        // value changed in return
            const label faci
        );

public:

    TypeName("CentredFitSGradData");


    // Constructors

        //- Construct from components
        CentredFitSGradData
        (
            const fvMesh& mesh,
            const extendedCentredCellToFaceStencil& stencil,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~CentredFitSGradData()
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
#   include "CentredFitSGradData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
