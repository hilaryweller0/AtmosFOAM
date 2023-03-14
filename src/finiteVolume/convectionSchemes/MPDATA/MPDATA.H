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
    Foam::fv::MPDATA

Description
    MPDATA is MPDATA with variable off-centred time stepping.
    Parameters set in a dictionary:
        nCorr 2;    // Number of correctors (default 1)
        gauge 1;    // MPDATA gauge mimimum (default 0)
        offCentre -1;   // "-1" implies set by Courant number or value [0,1]
        timeCorrector advective; // flux, advective or none
        FCTlimit false; // Apply Zalesak FCT, true or false
        FCTmin   0;     // If true give min and max or equal for monotone
        FCTmax   0;     // If true give min and max or equal for monotone
    
SourceFiles
    MPDATA.C

\*---------------------------------------------------------------------------*/

#ifndef MPDATA_H
#define MPDATA_H

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
                           Class MPDATA Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MPDATA
:
    public fv::convectionScheme<Type>
{
    // Private Data

        //- Dictionary for reading the MPDATA parameters
        dictionary dict_;

        //- Number of corrections (1 is normal)
        const label nCorr_;
        
        // MPDATA gauge (the minimum allowed value)
        const scalar gauge_;
        
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
        
        //- upwind interpolation using anti-diffusive flux
        tmp<gaussConvectionScheme<Type>> tantiDConvect_;

public:

    //- Runtime type information
    TypeName("MPDATA");

    // Constructors

        //- Construct from components
        MPDATA
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            label nCorr__ = 1,
            scalar gauge__ = 0,
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
            gauge_(gauge__),
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
            ),
            tantiDConvect_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        antiD(),
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, antiD())
                        )
                    )
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        MPDATA
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            dict_(is),
            nCorr_(readLabel(dict_.lookup("nCorr"))),
            gauge_(readScalar(dict_.lookup("gauge"))),
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
            ),
            tantiDConvect_
            (
                tmp<gaussConvectionScheme<Type>>
                (
                    new gaussConvectionScheme<Type>
                    (
                        mesh,
                        antiD(),
                        tmp<surfaceInterpolationScheme<Type>>
                        (
                            new upwind<Type>(mesh, antiD())
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
             && timeCorrector_ != "none"
            )
            {
                FatalErrorIn("MPDATA::MPDATA")
                    << "timeCorrector must be one of advective, flux or none but "
                    << timeCorrector_ << " was given" << exit(FatalError);
                
            }
        }
        
        //- Disallow default bitwise copy construction
        MPDATA(const MPDATA&) = delete;

        //- Destructor
        ~MPDATA()
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

        //- Access to anti-diffusive convection scheme
        const gaussConvectionScheme<Type>& antiDConvect() const
        {
            return tantiDConvect_();
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
        void operator=(const MPDATA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MPDATA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
