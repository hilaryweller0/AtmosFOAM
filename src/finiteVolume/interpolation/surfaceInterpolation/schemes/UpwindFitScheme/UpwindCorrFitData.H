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
    Foam::UpwindCorrFitData

Description
    Data for the quadratic fit correction interpolation scheme to be used with
    upwind biased stencil.
    - linearCorrection = true : fit calculated for corrected linear scheme
    - linearCorrection = false : fit calculated for corrected upwind scheme

SourceFiles
    UpwindCorrFitData.C

\*---------------------------------------------------------------------------*/

#ifndef UpwindCorrFitData_H
#define UpwindCorrFitData_H

#include "FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedUpwindCellToFaceStencilNew;

/*---------------------------------------------------------------------------*\
                    Class UpwindCorrFitData Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class UpwindCorrFitData
:
    public FitData
    <
        UpwindCorrFitData<Polynomial>,
        extendedUpwindCellToFaceStencilNew,
        Polynomial
    >
{
    // Private data

        //- For each face of the mesh store the coefficients to multiply the
        //  stencil cell values by if the flow is from the owner
        List<scalarList> owncoeffs_;

        //- For each face of the mesh store the coefficients to multiply the
        //  stencil cell values by if the flow is from the neighbour
        List<scalarList> neicoeffs_;


    // Private Member Functions

        //- Calculate the fit for the all the mesh faces
        //  and set the coefficients
        void calcFit();


public:

    TypeName("UpwindCorrFitData");


    // Constructors

        //- Construct from components
        UpwindCorrFitData
        (
            const fvMesh& mesh,
            const extendedUpwindCellToFaceStencilNew& stencil,
            const bool linearCorrection,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~UpwindCorrFitData()
    {}


    // Member functions

    //- Return reference to owner fit coefficients
    const List<scalarList>& owncoeffs() const
    {
        return owncoeffs_;
    }

    //- Return reference to neighbour fit coefficients
    const List<scalarList>& neicoeffs() const
    {
        return neicoeffs_;
    }

private:
    void fit(
            List<scalarList>& coeffs,
            const List<List<point> > stencilPoints,
            volScalarField& weightField,
            const label debugFaceI
    );
    void fit(
        const label faceI,
        List<scalarList>& coeffs,
        const List<List<point> > stencilPoints,
        const scalar wLin
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "UpwindCorrFitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
