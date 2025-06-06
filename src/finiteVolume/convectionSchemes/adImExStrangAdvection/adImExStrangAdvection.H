/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::fv::adImExStrangAdvection

Description
    Two time level adaptively implcit advection schemes with Strang
    carry-over splitting between off-centred trapezoidal implicit and an
    explicit RK scheme (typically RK3). 
    Parameters:
    alpha: variable off-centering
    beta : variable implicitness
    gamma: variable limiter for high-order flux correction
    are either constant or depend on the local Courant number.
    Schemes are a correction on 1st-order upwind.
    Parameters set in a dictionary:
        correctionScheme cubicUpwind; // Name of HO correction scheme
        RK_ButcherCoeffs 3((1 0 0) (0.25 0.250) (0.16667 0.16667 0.6666));
        alpha, beta and gamma (-1 for auto set based on Courant number)
SourceFiles
    adImExStrangAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef adImExStrangAdvection_H
#define adImExStrangAdvection_H

#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "dictionary.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class adImExStrangAdvection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class adImExStrangAdvection
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Dictionary for reading the adImExStrangAdvection parameters
        dictionary dict_;

        // Butcher tableau for explicit scheme as a square matrix
        const scalarSquareMatrix RK_;
        
        // Off centering uniform value [0,1] or variable (value <0)
        //const scalar alpha_;
        
        // Implicitness uniform value [0,1] or variable (value <0)
        //const scalar beta_;
        
        // HO limiter uniform value [0,1] or variable (value <0)
        //const scalar gamma_;
        
    // Private member functions
    
        // The upwind scheme
        tmp<surfaceInterpolationScheme<Type>> upwindScheme
        (
            const surfaceScalarField& advFlux
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                new upwind<Type>(this->mesh(), advFlux)
            );
        };
        // Upwind interpolation
        tmp<SurfaceField<Type>> upInterp
        (
            const surfaceScalarField& advFlux,
            const VolField<Type>& vf
        ) const
        {
            return upwindScheme(advFlux)->interpolate(vf);
        }
    
        // The correction scheme
        tmp<surfaceInterpolationScheme<Type>> corrScheme
        (
            const surfaceScalarField& advFlux
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                surfaceInterpolationScheme<Type>::New
                (
                    this->mesh(),
                    advFlux,
                    dict_.lookup("correctionScheme")
                )
            );
        };
        
        // The high-order correction
        tmp<SurfaceField<Type>> hCorr
        (
            const surfaceScalarField& advFlux,
            const VolField<Type>& vf
        ) const
        {
            return corrScheme(advFlux)->correction(vf);
        }

        
public:

    //- Runtime type information
    TypeName("adImExStrangAdvection");

    // Constructors

        //- Construct from mesh, advectingFlux and Istream
        adImExStrangAdvection
        (
            const fvMesh& mesh,
            const surfaceScalarField& advFlux,
            Istream& is
        );
        
        //- Disallow default bitwise copy construction
        adImExStrangAdvection(const adImExStrangAdvection&) = delete;

        //- Destructor
        ~adImExStrangAdvection() {}

    // Member Functions
    
        tmp<SurfaceField<Type>> interpolate
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        tmp<SurfaceField<Type>> flux
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        // The part of the advection which is solved implicitly
        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        tmp<VolField<Type>> fvcDiv
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adImExStrangAdvection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adImExStrangAdvection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
