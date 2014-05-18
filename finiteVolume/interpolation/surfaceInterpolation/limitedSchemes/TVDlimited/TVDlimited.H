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
    Foam::TVDlimited

Description
    TVDlimited limits the specified scheme to be TVD based on r obtained from
    the TVDNVD class.

SourceFiles
    TVDlimited.C

\*---------------------------------------------------------------------------*/

#ifndef TVDlimited_H
#define TVDlimited_H

#include "surfaceInterpolationScheme.H"
#include "limitedSurfaceInterpolationScheme.H"
#include "NVDTVD.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TVDlimited Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TVDlimited
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Reference to the surface flux used to choose upwind direction
        const surfaceScalarField& faceFlux_;

    // Private Member Functions

        //- Interpolation scheme
        tmp<surfaceInterpolationScheme<Type> > tInterp_;
        
        //- Disallow default bitwise copy construct
        TVDlimited(const TVDlimited&);

        //- Disallow default bitwise assignment
        void operator=(const TVDlimited&);


public:

    //- Runtime type information
    TypeName("TVDlimited");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        TVDlimited
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(mesh.lookupObject<surfaceScalarField>(word(is))),
            tInterp_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        TVDlimited
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            tInterp_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        inline tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        bool corrected() const { return false; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TVDlimitedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
