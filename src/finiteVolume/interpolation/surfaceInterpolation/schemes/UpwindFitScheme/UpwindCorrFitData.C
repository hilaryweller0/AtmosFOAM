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

\*---------------------------------------------------------------------------*/

#include "UpwindCorrFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "extendedUpwindCellToFaceStencilNew.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::UpwindCorrFitData<Polynomial>::UpwindCorrFitData
(
    const fvMesh& mesh,
    const extendedUpwindCellToFaceStencilNew& stencil,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        UpwindCorrFitData<Polynomial>,
        extendedUpwindCellToFaceStencilNew,
        Polynomial
    >
    (
        mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
    ),
    owncoeffs_(mesh.nFaces()),
    neicoeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing UpwindCorrFitData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<< "UpwindCorrFitData<Polynomial>::UpwindCorrFitData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::UpwindCorrFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    const label debugFaceI = 4388;
    volScalarField ownerWeights
    (
        IOobject
        (
            "ownerWeights",
       	    mesh.time().constant(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
        ),
	mesh,
	0.0,
	"fixedValue"
    );

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    // Owner stencil weights
    // ~~~~~~~~~~~~~~~~~~~~~

    // Get the cell/face centres in stencil order.
    List<List<point> > stencilPoints(mesh.nFaces());
    this->stencil().collectData
    (
        this->stencil().ownMap(),
        this->stencil().ownStencil(),
        mesh.C(),
        stencilPoints
    );

    // find the fit coefficients for every owner

    //Pout<< "-- Owner --" << endl;
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        scalarList wts(stencilPoints[facei].size(), scalar(1));
        FitData
        <
            UpwindCorrFitData<Polynomial>,
            extendedUpwindCellToFaceStencilNew,
            Polynomial
        >::calcFit(owncoeffs_[facei], wts, stencilPoints[facei], w[facei], facei);

        if (facei == debugFaceI)
        {
            forAll(stencilPoints[facei], stencilI)
            {
                forAll(mesh.C(), cellI)
                {
                    if (mesh.C()[cellI] == stencilPoints[facei][stencilI])
                    {
                        ownerWeights[cellI] = wts[stencilI];
                    }
                }
            }
        }

        //Pout<< "    facei:" << facei
        //    << " at:" << mesh.faceCentres()[facei] << endl;
        //forAll(owncoeffs_[facei], i)
        //{
        //    Pout<< "    point:" << stencilPoints[facei][i]
        //        << "\tweight:" << owncoeffs_[facei][i]
        //        << endl;
        //}
    }

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                scalarList wts(stencilPoints[facei].size(), scalar(1));
                FitData
                <
                    UpwindCorrFitData<Polynomial>,
                    extendedUpwindCellToFaceStencilNew,
                    Polynomial
                >::calcFit
                (
                    owncoeffs_[facei], wts, stencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }

    ownerWeights.write();


    // Neighbour stencil weights
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note:reuse stencilPoints since is major storage
    this->stencil().collectData
    (
        this->stencil().neiMap(),
        this->stencil().neiStencil(),
        mesh.C(),
        stencilPoints
    );

    // find the fit coefficients for every neighbour

    //Pout<< "-- Neighbour --" << endl;
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        scalarList wts(stencilPoints[facei].size(), scalar(1));
        FitData
        <
            UpwindCorrFitData<Polynomial>,
            extendedUpwindCellToFaceStencilNew,
            Polynomial
        >::calcFit(neicoeffs_[facei], wts, stencilPoints[facei], w[facei], facei);

        //Pout<< "    facei:" << facei
        //    << " at:" << mesh.faceCentres()[facei] << endl;
        //forAll(neicoeffs_[facei], i)
        //{
        //    Pout<< "    point:" << stencilPoints[facei][i]
        //        << "\tweight:" << neicoeffs_[facei][i]
        //        << endl;
        //}
    }

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                scalarList wts(stencilPoints[facei].size(), scalar(1));
                FitData
                <
                    UpwindCorrFitData<Polynomial>,
                    extendedUpwindCellToFaceStencilNew,
                    Polynomial
                >::calcFit
                (
                    neicoeffs_[facei], wts, stencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //
