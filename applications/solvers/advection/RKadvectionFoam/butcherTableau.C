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
    nSteps_(weights__.size()),
    subTimes_(subTimes__)
{
    if (nSteps_ != coeffs_.n())
    {
        FatalErrorIn("Foam::butcherTableau::butcherTableau from components")
            << " not all components have the same size:\n"
            << "coeffs_.n() = " << coeffs_.n() << nl
            << "weights_.size() = " << weights_.size() << nl
            << abort(FatalError);
    }
}


Foam::butcherTableau::butcherTableau(Istream& is)
:
    coeffs_(is),
    weights_(is),
    nSteps_(weights_.size()),
    subTimes_(is)
{
    if (nSteps_ != coeffs_.n())
    {
        FatalErrorIn("Foam::butcherTableau::butcherTableau from components")
            << " not all components have the same size:\n"
            << "coeffs_.n() = " << coeffs_.n() << nl
            << "weights_.size() = " << weights_.size()
            << abort(FatalError);
    }
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
