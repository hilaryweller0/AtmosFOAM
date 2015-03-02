/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    CLUST

Description
    CLUST interpolation scheme class derived from midPoint and
    returns blended linear/upwind weighting factors depending on the
    component of the wind through the face.

SourceFiles
    CLUST.C

\*---------------------------------------------------------------------------*/

#ifndef CLUST_H
#define CLUST_H

#include "midPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CLUST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CLUST
:
    public midPoint<Type>
{
    // Private Data

        //- The blending factor between linear (1-b) and linearUpwind (b)
        const scalar b_;

        // The vector field which defines the direction
        const surfaceVectorField& Uf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CLUST(const CLUST&);

        //- Disallow default bitwise assignment
        void operator=(const CLUST&);


public:

    //- Runtime type information
    TypeName("CLUST");


    // Constructors

        //- Construct from mesh and Istream
        CLUST
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
            b_(readScalar(schemeData)),
            Uf_
            (
                mesh.lookupObject<surfaceVectorField>
                (
                    word(schemeData)
                )
            )
        {
            if (!schemeData.eof())
            {
                IOWarningIn("CLUST(const fvMesh&, Istream&)", schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the face velocity field"
                       " should be specified."
                    << endl;
            }
        }

        //- Construct from mesh, faceFlux and Istream
        CLUST
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
            b_(readScalar(schemeData)),
            Uf_
            (
                mesh.lookupObject<surfaceVectorField>
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
