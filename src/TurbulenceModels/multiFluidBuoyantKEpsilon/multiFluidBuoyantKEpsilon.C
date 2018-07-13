/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "multiFluidBuoyantKEpsilon.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "Partitioned.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
multiFluidBuoyantKEpsilon<BasicTurbulenceModel>::multiFluidBuoyantKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool multiFluidBuoyantKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
multiFluidBuoyantKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    return
        (Cg_*this->Cmu_)*this->alpha_*this->k_*(g & fvc::grad(this->rho_))
       /(this->epsilon_ + this->epsilonMin_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
multiFluidBuoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    // kSource tmp to be returned
    tmp<fvScalarMatrix> kSourceTmp
    (
        new fvScalarMatrix(this->k_, dimensionSet(1,2,-3,0,0))
    );
    fvScalarMatrix& kS = kSourceTmp.ref();

    const fvMesh& mesh = this->mesh_;
    
    const uniformDimensionedVectorField& g =
        mesh.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    const 
     Partitioned<ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>>&
       turbulence =
       Partitioned<ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>>::New
       (mesh);

    const word& kName = this->k_.name();
    const wordList& partNames = turbulence.partNames();

    // Get the name and number of this partition
    const word partName = kName.substr(kName.find(".") + 1);
    const label ip = findIndex(partNames, partName);

    // Accumulate the diaganol and the RHS of the matrix for the source term
    volScalarField diag
    (
        IOobject("kSourceDiag", mesh.time().timeName(), mesh), mesh,
        dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), scalar(0))
    );
    // Loop through the different partitions
    for(label jp = 0; jp < turbulence.size(); jp++)
    {
        // Other partition name
        const word& otherPartName = partNames[jp];
    
        // k for the other partition
        const volScalarField& kj = mesh.objectRegistry::template
              lookupObject<volScalarField>("k."+partNames[jp]);
    
        // Check that it is not this partition
        if (ip != jp)
        {
            const volScalarField& massTransferij
                 = mesh.objectRegistry::template
      lookupObject<volScalarField>("massTransfer."+partName+'.'+partNames[jp]);
            const volScalarField& massTransferji
                 = mesh.objectRegistry::template
      lookupObject<volScalarField>("massTransfer."+partNames[jp]+'.'+partName);
        
            diag += massTransferij;
            
            kS.source() += massTransferji*kj;
        }
    }

    if (mag(g.value()) > small)
    {
        kS = -fvm::SuSp(diag+Gcoef(), this->k_);
    }
    else
    {
        kS = - fvm::SuSp(diag, this->k_)
             + kEpsilon<BasicTurbulenceModel>::kSource();
    }
    return kSourceTmp;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
multiFluidBuoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        vector gHat(g.value()/mag(g.value()));

        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("small", dimVelocity, small)
        );

        return -fvm::SuSp(this->C1_*tanh(mag(v)/u)*Gcoef(), this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::epsilonSource();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
