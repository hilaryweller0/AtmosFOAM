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
    Foam::fv::adaptiveImplicitAdvection

Description
    Two time level advection schemes with variable off-centering, alpha, and
    variable implicitness depending on CoLimit.
    Schemes are a correction on 1st-order upwind.
    Parameters set in a dictionary:
        correctionScheme QUICKupwind 1 2; // Name of HO correction scheme
        nCorr 1;    // Number of correctors (default 1)
        offCentre [0,1] or <0; // <0 means calculated from Courant no, Co
        CoLimit   1; // Max Co to be treated explicitly
        //fullSolver true; // Use full solver (true) or HbyA
        FCTlimit false; // Apply Zalesak FCT, true or false
        FCTmin   0;     // If true give min and max or equal for monotone
        FCTmax   0;     // If true give min and max or equal for monotone
    The faceFlux is assumed to be at the middle of the time level
    
SourceFiles
    adaptiveImplicitAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveImplicitAdvection_H
#define adaptiveImplicitAdvection_H

#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveImplicitAdvection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class adaptiveImplicitAdvection
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Dictionary for reading the adaptiveImplicitAdvection parameters
        dictionary dict_;

        //- Number of HO corrections
        const int nCorr_;

        // Off centering uniform value [0,1] or variable (value <0)
        const scalar offCentre_;
        
        // Max Courant number to be treated explicitly
        const scalar CoLimit_;
        
        // Full solver or HbyA solution
        //const Switch fullSolver_;
        
        // Apply Zalesak FCT correction with given limits (equal for monotone)
        const Switch FCTlimit_;
        const scalar FCTmin_, FCTmax_;

    // Private member functions
    
        //- Calculate the off centering that varies in space
        tmp<surfaceScalarField> calcOffCentre
        (
            const surfaceScalarField& flux
        ) const;
        
        // The correction scheme
        tmp<surfaceInterpolationScheme<Type>> corrScheme
        (
            const surfaceScalarField& flux
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                surfaceInterpolationScheme<Type>::New
                (
                    this->mesh(),
                    flux,
                    dict_.lookup("correctionScheme")
                )
            );
        };

        
public:

    //- Runtime type information
    TypeName("adaptiveImplicitAdvection");

    // Constructors

        //- Construct from mesh, faceFlux and Istream
        adaptiveImplicitAdvection
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );
        
        //- Disallow default bitwise copy construction
        adaptiveImplicitAdvection(const adaptiveImplicitAdvection&) = delete;

        //- Destructor
        ~adaptiveImplicitAdvection() {}

    // Member Functions
    
        // Return the offCentre field
        const surfaceScalarField& offCentre() const ;
        
        // Maximum Courant number to be treated explicitly
        scalar CoLimit() const { return CoLimit_; }
    
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        // The part of the advection which is solved implicitly
        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adaptiveImplicitAdvection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adaptiveImplicitAdvection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
