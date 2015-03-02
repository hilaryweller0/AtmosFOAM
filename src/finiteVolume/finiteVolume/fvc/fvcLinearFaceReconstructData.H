/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fvcLinearFaceReconstructData

Description
    Data for the linear fit for reconstructing face vectors from face fluxes

SourceFiles
    fvcLinearFaceReconstructData.C

\*---------------------------------------------------------------------------*/

#ifndef fvcLinearFaceReconstructData_H
#define fvcLinearFaceReconstructData_H

#include "MeshObject.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class centredCFCFaceToFaceStencilObject;

/*---------------------------------------------------------------------------*\
                    Class fvcLinearFaceReconstructData Declaration
\*---------------------------------------------------------------------------*/

class fvcLinearFaceReconstructData
:
    public MeshObject<fvMesh, fvcLinearFaceReconstructData>
{
    // Private data

        //- The stencil the fit is based on
        const centredCFCFaceToFaceStencilObject& stencil_;

        //- Dimensionality of the geometry
        const label dim_;

        //- For spherical coords, the k direction should be outwards
        const bool spherical_;

        //- The number of valid (ie non-empty) faces, internal and boundary
        const label nValidFaces_;

        //- For each face in the mesh store the values which multiply the
        //  values of the stencil to obtain the reconstructed vector
        List<List<vector> > coeffs_;


    // Private member functions

        //- Add the coefficients to a line of the matrix to be inverted based on
        //  the values for one point of the stencil
        void addCoeffs
        (
            scalar* coeffs,
            const vector& d,
            const vector& Sf,
            const label faceI,
            const bool smallStencil
        );

        //- Find the i, j and k directions for each face
        void findFaceDirs
        (
            vector& idir,         // value changed in return
            vector& jdir,         // value changed in return
            vector& kdir,         // value changed in return
            const label faceI
        );


public:

    // TypeName("fvcLinearFaceReconstructData");

    // Constructors

        //- Construct from components
        fvcLinearFaceReconstructData
        (
            const fvMesh& mesh,
            const centredCFCFaceToFaceStencilObject& stencil
        );


    //- Destructor
    virtual ~fvcLinearFaceReconstructData()
    {}


    // Member functions

        //- Return reference to fit coefficients
        const List<List<vector> >& coeffs() const
        {
            return coeffs_;
        }

        //- Return reference to the stencil
        const centredCFCFaceToFaceStencilObject& stencil() const
        {
            return stencil_;
        }

        //- Calculate the fit for the specified cell and set the coefficients
        void calcFit
        (
            List<vector>& coeffsi, // coefficients to be set
            const List<point>&,    // Stencil points
            const List<vector>&,   // Stencil face area vectors
            const label faceI      // Current face index
        );

        //- Calculate the fit for all the cells
        void calcFit();

        //- Delete the data when the mesh moves not implemented
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvcLinearFaceReconstructData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
