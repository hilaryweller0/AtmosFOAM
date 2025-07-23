/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2025 OpenFOAM Foundation
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
    Foam::functionObjects::multiScalarRKTransport

Description
    Transports a set of scalars, possibly including density

    Example:
    \verbatim
        #includeFunc multiScalarRKTransport
        (
            nFields=nFields,
            fields=(rho T),
            RK_ButcherCoeffs=3 3((1 0 0)
                             (0.25 0.250 0)
                             (0.16667 0.16667 0.66666))
        )
    \endverbatim

    with scheme specification:
        div(phiv,rho)     Gauss linear;
        multiScalarRKTransport
        {
            fields (rho T);
            advectingFlux phiv;
            massFlux [none|phi];
            RK_ButcherCoeffs 3 3((1 0 0)
                             (0.25 0.250 0)
                             (0.16667 0.16667 0.66666));
            gammaCoeffs (2.5 6.5);
            correctionScheme cubicUpwind;
        }

    and solver specification:
        multiScalarRKTransport
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0;
        }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    multiScalarRKTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_multiScalarRKTransport_H
#define functionObjects_multiScalarRKTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "scalarMatrices.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class multiScalarRKTransport Declaration
\*---------------------------------------------------------------------------*/

class multiScalarRKTransport
:
    public fvMeshFunctionObject
{

private:

    // Private Data

        //- Number of fields being advected
        const int nFields_;

        //- Names of fields to process
        wordList fieldNames_;

        //- The fields being advected
        PtrList<volScalarField> s_;

        //- Name of advecting flux
        word phivName_;

        //- Name of calculated mass flux field (none if no density)
        word massFluxName_;

        //- Butcher tableau for explicit scheme as a square matrix
        const scalarSquareMatrix RK_;

        //- Coefficients for calculating gamma from the Courant number
        const FixedList<scalar,2> gammaCoeffs_;

        //- Name of the high-order correction scheme
        const word correctionSchemeName_;

        //- How many FCT iterations to apply for each field (0 if no FCT)
        const labelList FCTiter_;

        //- (Optional) min and max values for each tracer
        scalarListList FCTlimits_;
        
        //- The dictionary that describes the (changing) velocity field
        const word velocityDictName_;
        
        //- The advecting flux
        surfaceScalarField phiv_;

    // Private member functions

        // The upwind scheme
        tmp<surfaceInterpolationScheme<scalar>> upwindScheme
        (
            const surfaceScalarField& flux
        ) const
        {
            return tmp<surfaceInterpolationScheme<scalar>>
            (
                new upwind<scalar>(mesh_, flux)
            );
        };
        // Upwind interpolation
        tmp<SurfaceField<scalar>> upInterp
        (
            const surfaceScalarField& flux,
            const volScalarField& vf
        ) const
        {
            return upwindScheme(flux)->interpolate(vf);
        }

        // The correction scheme
        tmp<surfaceInterpolationScheme<scalar>> corrScheme
        (
            const surfaceScalarField& flux
        ) const
        {
            IStringStream is(correctionSchemeName_);
            return tmp<surfaceInterpolationScheme<scalar>>
            (
                surfaceInterpolationScheme<scalar>::New(mesh_, flux, is)
            );
        };

        // The high-order correction
        tmp<SurfaceField<scalar>> hCorr
        (
            const surfaceScalarField& flux,
            const volScalarField& vf
        ) const
        {
            return corrScheme(flux)->correction(vf);
        }
        
        // Set the flux for the wind at an intermediate time. i=-1 for old or
        // new timelevel, i>=0 for an RK stage
        void setFlux(const label i=-1);

public:

    //- Runtime type information
    TypeName("multiScalarRKTransport");


    // Constructors

        //- Construct from Time and dictionary
        multiScalarRKTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multiScalarRKTransport(const multiScalarRKTransport&) = delete;


    //- Destructor
    virtual ~multiScalarRKTransport();


    // Member Functions

        //- Read the multiScalarRKTransport data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do not execute at the start of the run
        virtual bool executeAtStart() const
        {
            return false;
        }

        //- Calculate the multiScalarRKTransport
        virtual bool execute();

        //- Write the updated scalar field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiScalarRKTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
