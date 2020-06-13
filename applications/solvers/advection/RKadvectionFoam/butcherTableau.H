/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::butcherTableau

Description
    Holds the coefficients and weights for a Butcher Tableau
    for an explicit Runge-Kutta scheme. Member functions sum fields for
    RK steps and for the final step.

SourceFiles
    butcherTableau.C

\*---------------------------------------------------------------------------*/

#ifndef butcherTableau_H
#define butcherTableau_H

#include "scalarMatrices.H"
#include "RKfield.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                         Class butcherTableau Declaration
\*---------------------------------------------------------------------------*/

class butcherTableau
{
    // Private data

        //- RK coefficients
        const scalarSquareMatrix coeffs_;
        
        //- RK weights (set to zero if final row to coeffs to be used instead)
        const scalarList weights_;

        //- Number of RK steps
        const label nSteps_;
        
        //- Sub-times (times for each row of coeffs)
        const scalarList subTimes_;

    // Private Member Functions

public:

    // Constructors

        //- Construct from components
        butcherTableau
        (
            const scalarSquareMatrix& coeffs__,
            const scalarList& weights__,
            const scalarList& subTimes__=scalarList(0)
        );

        //- Construct from Istream
        butcherTableau(Istream&);

        //- Construct as copy
        butcherTableau(const butcherTableau&);


    // Member Functions

        // Access
        const scalarSquareMatrix& coeffs() const { return coeffs_; }
        const scalarList& weights() const { return weights_; }
        label nSteps() const { return nSteps_; }
        const scalarList& subTimes() const { return subTimes_; }

        // Use
        
        //- Sum fields for a single RK step
        template<class Type,template<class> class PatchField,class GeoMesh>
        tmp<GeometricField<Type, PatchField, GeoMesh> >
        RKstep
        (
            const label i,
            const RKfield<Type, PatchField, GeoMesh>& rkField,
            label jmax = -1
        ) const;

        //- Sum fields for the final RK step
        template<class Type,template<class> class PatchField,class GeoMesh>
        tmp<GeometricField<Type, PatchField, GeoMesh> >
        RKfinal
        (
            const RKfield<Type, PatchField, GeoMesh>& rkField
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "butcherTableau.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
