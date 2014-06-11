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
    Foam::LUSTbb

Description
    LUSTbb: Linear-upwind stabilised transport with a blending factor, b

    Interpolation scheme class derived from linearUpwind which returns blended
    linear/linear-upwind weighting factors and also applies a explicit
    gradient-based correction obtained from the linearUpwind scheme.  The
    blending-factor is set to b linear which optimises the balance between
    accuracy and stability on a range of LES cases with a range of mesh quality.

SourceFiles
    LUSTb.C

\*---------------------------------------------------------------------------*/

#ifndef LUSTb_H
#define LUSTb_H

#include "midPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LUSTb Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LUSTb
:
    public midPoint<Type>
{
    // Private Data
    
        //- The blending factor between linear (1-b) and linearUpwind (b)
        const scalar b_;

        // the flux field to set the upwind directin
        const surfaceScalarField& un_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LUSTb(const LUSTb&);

        //- Disallow default bitwise assignment
        void operator=(const LUSTb&);


public:

    //- Runtime type information
    TypeName("LUSTb");


    // Constructors

        //- Construct from mesh and Istream
        LUSTb
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
            b_(readScalar(schemeData)),
            un_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(schemeData)
                )
            )
        {
            if (!schemeData.eof())
            {
                IOWarningIn("LUSTb(const fvMesh&, Istream&)", schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the normal velocity "
                       " should be specified."
                    << endl;
            }
        }

        //- Construct from mesh, faceFlux and Istream
        LUSTb
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
            b_(readScalar(schemeData)),
            un_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(schemeData)
                )
            )
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
