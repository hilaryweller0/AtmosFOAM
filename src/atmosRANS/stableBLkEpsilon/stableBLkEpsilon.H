/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::RASModels::stableBLkEpsilon

Description
    Additional buoyancy generation/dissipation term applied to the
    k and epsilon equations of the standard k-epsilon model.

    References:
    \verbatim
        Duynkerke P.G. (1988).
        Application of the {E--$\varepsilon$} Turbulence Closure Model to the
        Neutral and Stable Atmospheric Boundary Layer.
        J. Atmos. Sci. 45(5):865-880
        
        Cuxart et al (2006).
        Single-Column Model Intercomparison for a Stably Stratified Atmospheric
        Boundary Layer.
        Boundary Layer Meteorology. 118(2):273-303
    \endverbatim

    Modifications based on testing.

    The default model coefficients are
    \verbatim
        stableBLkEpsilon
        {
            Cg              1.0;
            kappa           0.41;
            betam           4.8;
            betah           7.8;
            Tref            288;
        }
    \endverbatim

SourceFiles
    stableBLkEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef stableBLkEpsilon_H
#define stableBLkEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class stableBLkEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class stableBLkEpsilon
:
    public kEpsilon<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cg_;
            dimensionedScalar kappa_;
            dimensionedScalar betam_;
            dimensionedScalar betah_;
            dimensionedScalar Tref_;

        // Surface fields
            volScalarField z0_;
            volScalarField LmoInv_;
            volScalarField uStar_;
            volScalarField thetaStar_;
            volScalarField nutGround_;
            volScalarField alphatGround_;
            volScalarField epsilonGroundCorr_;
            
        // Vol fields
            volScalarField theta_;
        
        // Themodynamic pRef and R/Cp for updating potential temperature
            const dimensionedScalar& pRef_;
            scalar RbyCp_;

    // Protected Member Functions

        //- Bound epsilon and return Cmu*sqr(k) for nut
        tmp<volScalarField> boundEpsilon();

        //- Correct the eddy-viscosity nut
        virtual void correctNut();

        //- Source term for the k equation
        virtual tmp<fvScalarMatrix> kSource() const;

        //- Source term for the epsilon equation
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        // Bouyancy source
        tmp<volScalarField> Gcoef() const;
        
        // Update the surface Fields
        virtual void updateSurfaceFields();

        // Update theta
        volScalarField& calcTheta();

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("stableBLkEpsilon");


    // Constructors

        //- Construct from components
        stableBLkEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        stableBLkEpsilon(const stableBLkEpsilon&) = delete;


    //- Destructor
    virtual ~stableBLkEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        //- and update the surface fields
        virtual void correct();

    // Access
    
        const dimensionedScalar& Cmu() const { return this->Cmu_; }
        const dimensionedScalar& kappa() const { return this->kappa_; }
    
        //- The surface roughness
        const volScalarField& z0() const { return z0_;}

        //- The Monin-Obukov length
        const volScalarField& LmoInv() const { return LmoInv_;}

        //- The friction velocity
        const volScalarField& uStar() const { return uStar_;}

        //- The surface thetaStar
        const volScalarField& thetaStar() const { return thetaStar_;}

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stableBLkEpsilon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stableBLkEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
