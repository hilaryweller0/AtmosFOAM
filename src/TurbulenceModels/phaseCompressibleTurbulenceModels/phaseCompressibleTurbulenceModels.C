/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "PhaseCompressibleTurbulenceModel.H"

#include "transportModel.H"
#include "turbulentFluidThermoModel.H"
#include "rhoThermo.H"

#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    fluidThermo
);

makeBaseTurbulenceModel
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    fluidThermo
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (fluidThermoPhaseCompressibleTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoPhaseCompressibleTurbulenceModel, LES, Type)

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "multiFluidBuoyantKEpsilon.H"
makeRASModel(multiFluidBuoyantKEpsilon);

#include "multiFluidThetaBuoyantKEpsilon.H"
makeRASModel(multiFluidThetaBuoyantKEpsilon);


// ************************************************************************* //
