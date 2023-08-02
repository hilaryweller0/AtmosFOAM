/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InNamespace
    Foam::fv

Description
    Applies Zalesak (1979) flux limiter to an advection-diffusion equation 
    and solves. 

SourceFiles
    fvFCTadvectionDiffusion.C

\*---------------------------------------------------------------------------*/


#ifndef fvFCTadvectionDiffusion_H
#define fvFCTadvectionDiffusion_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fv
{
    // Advect and diffuse T, density weighted with rho, offCentre varying in 
    // space, with flux, diffusivity gamma and explicit source term S
    void FCTadvectionDiffusion
    (
        volScalarField& T,
        const volScalarField& rho,
        const surfaceScalarField& offCentre,
        const surfaceScalarField& flux,
        const surfaceScalarField& gamma,
        const volScalarField& S,
        const bool implicitAdvection = true,
        const bool finalIter = true
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
