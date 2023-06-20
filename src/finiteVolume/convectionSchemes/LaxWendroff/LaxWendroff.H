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
    Foam::fv::LaxWendroff

Description
    LaxWendroff with variable off-centred time stepping as an explicit
     correction on first-order upwind (implicit or explicit)
    Parameters set in a dictionary:
        offCentre -1;   // "-1" implies set by Courant number or value [0,1]
        FCTlimit false; // Apply Zalesak FCT, true or false
        FCTmin   0;     // If true give min and max or equal for monotone
        FCTmax   0;     // If true give min and max or equal for monotone
    
SourceFiles
    LaxWendroff.C

\*---------------------------------------------------------------------------*/

#ifndef LaxWendroff_H
#define LaxWendroff_H

#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "dictionary.H"
#include "orthogonalSnGrad.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class LaxWendroff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LaxWendroff
:
    public fv::convectionScheme<Type>
{
    // Private Data

        // Dictionary for reading the LaxWendroff parameters
        dictionary dict_;

        // Off centering uniform value [0,1] or variable (value <0)
        const scalar offCentre_;
        
        //- Use advective or flux form time corrector or none
        word timeCorrector_;

        //- Apply Zalesak FCT correction with given limits (equal for monotone)
        const Switch FCTlimit_;
        const Type FCTmin_, FCTmax_;

        // Density field name
        const word densityName_;

        // Schemes needed 

        // upwind divergence
        gaussConvectionScheme<Type> GaussUpwind;
        
        // Max interpolate
        localMax<scalar> maxInterp;
        
        
    // Private member functions
    
        //- Calculate the off centering that varies in space
        tmp<surfaceScalarField> calcOffCentre
        (
            const surfaceScalarField& flux
        ) const;
        
        //- Smooth the correction where offCentre > 0
        void smoothFluxCorr
        (
            surfaceScalarField& fluxCorr,
            const surfaceScalarField& offCentre
        ) const;
        
        //- Calculate the high order flux correction
        void calculateCorr
        (
            surfaceScalarField& Tcorr,
            const surfaceScalarField& volFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& offCentre
        ) const;

        void advect
        (
            volScalarField& T,
            const surfaceScalarField& faceFlux,
            const surfaceScalarField& volFlux
        ) const;
                
        // Return the density field
        const volScalarField& rho() const
        {
            return this->mesh().objectRegistry::template
                     lookupObject<volScalarField>(densityName_);
        }
        
        // Check for density field
        bool isRho() const
        {
            return this->mesh().objectRegistry::template
                 foundObject<volScalarField>(densityName_);
        }
        
        // Divide by the mid-time density on the face
        tmp<surfaceScalarField> specific(const surfaceScalarField& flux) const;

        // Divide by the density
        tmp<GeometricField<Type, fvPatchField, volMesh>> 
        specific(const GeometricField<Type, fvPatchField, volMesh>& vf) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> 
        specific
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
        ) const;
        
        tmp<surfaceVectorField> upwindDelta
        (
            const surfaceScalarField& flux
        ) const;


public:

    //- Runtime type information
    TypeName("LaxWendroff");

    // Constructors

        //- Construct from mesh, faceFlux and Istream
        LaxWendroff
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );
        
        //- Disallow default bitwise copy construction
        LaxWendroff(const LaxWendroff&) = delete;

        //- Destructor
        ~LaxWendroff(){}

    // Member Functions
    
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
        void operator=(const LaxWendroff&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LaxWendroff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
