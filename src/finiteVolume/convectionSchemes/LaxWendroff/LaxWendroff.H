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
        nCorr 2;    // Number of correctors (default 1)
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

        //- Dictionary for reading the LaxWendroff parameters
        dictionary dict_;

        //- Number of corrections (1 is normal)
        const label nCorr_;
        
        // Off centering uniform value [0,1] or variable (value <0)
        const scalar offCentre_;
        
        // Use advective or flux form time corrector or none
        word timeCorrector_;

        // Apply Zalesak FCT correction with given limits (equal for monotone)
        const Switch FCTlimit_;
        const Type FCTmin_, FCTmax_;
        
        const surfaceScalarField faceVol_;

        //- anti-diffusive flux
        mutable surfaceScalarField* antiDPtr_;

        //- upwind interpolation
        tmp<gaussConvectionScheme<Type>> tupwindConvection_;

public:

    //- Runtime type information
    TypeName("LaxWendroff");

    // Constructors

        //- Construct from components
        LaxWendroff
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            label nCorr__ = 1,
            scalar offCentre__ = 0,
            word timeCorrector__ = "advective",
            Switch FCTlimit__ = false,
            scalar FCTmin__ = 0,
            scalar FCTmax__ = 0
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            dict_(),
            nCorr_(nCorr__),
            offCentre_(offCentre__),
            timeCorrector_(timeCorrector__),
            FCTlimit_(FCTlimit__),
            FCTmin_(FCTmin__),
            FCTmax_(FCTmax__),
            faceVol_(mesh.Sf() & mesh.delta()),
            antiDPtr_(new surfaceScalarField(faceFlux)),
            tupwindConvection_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        faceFlux,
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, faceFlux)
                        )
                    )
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        LaxWendroff
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            dict_(is),
            nCorr_(readLabel(dict_.lookup("nCorr"))),
            offCentre_(readScalar(dict_.lookup("offCentre"))),
            timeCorrector_(dict_.lookupOrDefault<word>("timeCorrector", "advective")),
            FCTlimit_(dict_.lookup("FCTlimit")),
            FCTmin_(dict_.lookupOrDefault<scalar>("FCTmin", scalar(0))),
            FCTmax_(dict_.lookupOrDefault<scalar>("FCTmax", scalar(0))),
            faceVol_(mesh.Sf() & mesh.delta()),
            antiDPtr_(new surfaceScalarField(faceFlux)),
            tupwindConvection_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        faceFlux,
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, faceFlux)
                        )
                    )
                )
            )
        {
            if (offCentre_ >= 0.5 - SMALL) timeCorrector_ = "none";
        
            else if
            (
                timeCorrector_ != "advective"
             && timeCorrector_ != "flux"
             && timeCorrector_ != "RK2"
             && timeCorrector_ != "none"
            )
            {
                FatalErrorIn("LaxWendroff::LaxWendroff")
                    << "timeCorrector must be one of advective, flux or none but "
                    << timeCorrector_ << " was given" << exit(FatalError);
                
            }
        }
        
        //- Disallow default bitwise copy construction
        LaxWendroff(const LaxWendroff&) = delete;

        //- Destructor
        ~LaxWendroff()
        {
            delete antiDPtr_;
        }

    // Member Functions
    
        // Return the anti-diffusive flux
        surfaceScalarField& antiD() const { return *antiDPtr_; }
    
        //- Access to upwind convection scheme
        const gaussConvectionScheme<Type>& upwindConvect() const
        {
            return tupwindConvection_();
        }

        //- Calculate the anti-diffusive flux
        void calculateAntiD
        (
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& offCentre
        ) const;

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
