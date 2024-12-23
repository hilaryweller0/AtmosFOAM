/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#ifndef fvCFD_H
#define fvCFD_H

#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "fvcDiv.H"
#include "linear.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "fvmLaplacian.H"
#include "fvmDiv.H"
#include "fvcAverage.H"
#include "fvcCellReduce.H"
#include "fvcCurl.H"
#include "fvcD2dt2.H"
#include "fvcDdt.H"
#include "fvcDDt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvc.H"
#include "fvcLaplacian.H"
#include "fvcMagSqrGradGrad.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "fvcSup.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcVolumeIntegrate.H"


#include "fvScalarMatrix.H"
#include "fvModels.H"
using namespace Foam;

#endif
// ************************************************************************* //
