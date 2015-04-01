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
    Foam::APVM

Description
    Generalised "Anticipated Potential Vorticity Method" for interpolation 
    of any volField onto a surfaceField. For a volField Psi, the
    interpolation is:
    midPoint(Psi) + 0.5*dt Uf . grad(Psi)
    where Uf is the face velocity, . is the dot product and grad(Psi)
    is the face gradient. Uf and grad are looked up from the database.
    The scheme returns an explicit correction on midPoint

SourceFiles
    APVM.C

\*---------------------------------------------------------------------------*/

#ifndef APVM_H
#define APVM_H

#include "midPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class APVM Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class APVM
:
    public midPoint<Type>
{
    // Private Data

        const surfaceVectorField& Uf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        APVM(const APVM&);

        //- Disallow default bitwise assignment
        void operator=(const APVM&);


public:

    //- Runtime type information
    TypeName("APVM");


    // Constructors

        //- Construct from Istream.
        //  The name of the velocity field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        APVM
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
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
                IOWarningIn("APVM(const fvMesh&, Istream&)", schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the face velocity field"
                       " should be specified."
                    << endl;
            }
        }

        //- Construct from mesh, surfaceScalarField and Istream.
        APVM
        (
            const fvMesh& mesh,
            const surfaceScalarField& Uf,
            Istream& schemeData
        )
        :
            midPoint<Type>(mesh),
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
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
