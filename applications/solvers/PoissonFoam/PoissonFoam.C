/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    PoissonFoam

Description
    Solve the elliptic Poisson equation, laplacian(T) = RHS

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const dictionary& solutionDict =
        mesh.solutionDict().subDict("potentialFlow");

    const int nCorr = solutionDict.lookupOrDefault<int>("nCorrectors", 0);
    const label TRefCell = solutionDict.lookupOrDefault<label>
    (
        "referenceCell", 0
    );
    const scalar TRefValue = solutionDict.lookupOrDefault<scalar>
    (
        "referenceValue", scalar(0)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Increment time so that results are written out at a different
    // time to the start time
    runTime++;
    
    // Outer iterations for a non-orthogonal mesh (so laplacian is not
    // discretised fully implicitly)
    for(int iCorr = 0; iCorr < nCorr; iCorr++)
    {
        // Setup matrix equation and solve
        fvScalarMatrix TEqn(fvm::laplacian(T) == RHS);
        TEqn.setReference(TRefCell, TRefValue);
        TEqn.solve();
    }

    runTime.write();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
