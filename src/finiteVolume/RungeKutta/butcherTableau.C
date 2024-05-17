/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "butcherTableau.H"
using namespace Foam;

// * * * * * * * * * * * * * Private Member Function * * * * * * * * * * * //

void Foam::butcherTableau::checkAndFix()
{
    // Enforce subTimes to be the sum of the rows of coeffs
    subTimes_.resize(nSteps_);
    for(label i = 0; i < nSteps_; i++)
    {
        subTimes_[i] = 0;
        for(label j = 0; j < nSteps_; j++) subTimes_[i] += coeffs_[i][j];
    }
    
    // Check that the weights are zero or the correct length sum to one
    if (weights_.size() > 0)
    {
        if (weights_.size() != nSteps_ && max(mag(weights_)) > SMALL)
        {
            FatalErrorIn("Foam::butcherTableau::checkAndFix")
            << " Butcher tableau is " << coeffs_.n() << " by " << coeffs_.n()
            << " weights are " << weights_
            << abort(FatalError);
        }
        else if (weights_.size() == nSteps_ && max(mag(weights_)) > SMALL)
        {
            // Check that the weights sum to one
            scalar wSum = 0;
            for(label i = 0; i < nSteps_; i++) wSum += weights_[i];
            if (mag(wSum - 1) > SMALL)
            {
                FatalErrorIn("Foam::butcherTableau::checkAndFix")
                    << " weights are " << weights_ << " which sum to " << wSum
                    << " but which should sum to 1" << abort(FatalError);
            }
        }
        else if (weights_.size() != nSteps_ && max(mag(weights_)) <= SMALL)
        {
            weights_.resize(nSteps_);
            weights_ = scalar(0);
        }
        //else  (weights_.size() == nSteps_ && max(mag(weights_)) <= SMALL)
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::butcherTableau::butcherTableau
(
    const scalarSquareMatrix& coeffs__,
    const scalarList& weights__,
    const scalarList& subTimes__
)
:
    coeffs_(coeffs__),
    weights_(weights__),
    nSteps_(coeffs__.n()),
    subTimes_(subTimes__)
{
    checkAndFix();
}


Foam::butcherTableau::butcherTableau(Istream& is)
:
    coeffs_(is),
    weights_(is),
    nSteps_(coeffs_.n()),
    subTimes_(is)
{
    checkAndFix();
}


Foam::butcherTableau::butcherTableau(const butcherTableau& bt)
:
    coeffs_(bt.coeffs_),
    weights_(bt.weights_),
    nSteps_(bt.nSteps_),
    subTimes_(bt.subTimes_)
{}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<GeometricField<Type, PatchField, GeoMesh> >
Foam::butcherTableau::RKstep
(
    const label i,
    const RKfield<Type, PatchField, GeoMesh>& rkField,
    label jmax
) const
{
    if (jmax == -1) jmax = i-1;
    
    tmp<GeometricField<Type, PatchField, GeoMesh> > tStep
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "RKstep",
                rkField[0].instance(),
                rkField[0].mesh()
            ),
            rkField[0].mesh(),
            dimensioned<Type>("rk",rkField[0].dimensions(),0*rkField[0][0])
        )
    );
    
    for(label j = 0; j <= jmax; j++)
    {
        tStep.ref() == tStep.ref() + coeffs_[i][j]*rkField[j];
    }
    
    return tStep;
}

template<class Type,template<class> class PatchField,class GeoMesh>
Foam::tmp<GeometricField<Type, PatchField, GeoMesh> >
Foam::butcherTableau::RKfinal
(
    const RKfield<Type, PatchField, GeoMesh>& rkField
) const
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tFinal
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "RKfinal",
                rkField[0].instance(),
                rkField[0].mesh()
            ),
            rkField[0].mesh(),
            dimensioned<Type>("rk",rkField[0].dimensions(),0*rkField[0][0])
        )
    );
    
    for(label i = 0; i < nSteps(); i++)
    {
        tFinal.ref() == tFinal.ref() + weights_[i]*rkField[i];
    }
    
    return tFinal;

}

// ************************************************************************* //
