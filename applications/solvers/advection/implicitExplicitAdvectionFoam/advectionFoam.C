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
    advectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "velocityField.H"
#include "CourantNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()

    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<label>("nCorr", label(2));
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    const scalar CoLimit = readScalar(mesh.schemesDict().lookup("CoLimit"));


    // Output file to write error measures each time step
    OFstream os("errorMeasures.dat");
    os << "#Time minT maxT sumT TV" << endl;

    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "advectionDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    bool timeVaryingWind = dict.lookupOrDefault<bool>("timeVaryingWind", false);
    const dictionary& velocityDict = dict.subOrEmptyDict("velocity");
    autoPtr<velocityField> v;
    if (velocityDict.size() > 0)
    {
        v = velocityField::New(dict.subOrEmptyDict("velocity"));
    }
    
    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        if (timeVaryingWind)
        {
            v->applyTo(phi);

            forAll(phi, faceI)
            {
                if (mag(phi[faceI]) > phiLimit[faceI])
                {
                    phiSmall[faceI] = 0;
                    phiBig[faceI] = phi[faceI];
                }
                else
                {
                    phiSmall[faceI] = phi[faceI];
                    phiBig[faceI] = 0;
                }
            }

            U = fvc::reconstruct(phi);
            Uf = linearInterpolate(U);
            Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        }

        #include "CourantNo.H"

        for(label corr = 0; corr < nCorr; corr++)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + (1-offCentre)*divPhiOld
              + offCentre*fvc::div(phiSmall, T, "explicit")
              + offCentre*fvm::div(phiBig, T, "upwind")
              - offCentre*fvc::div(phiBig, T, "upwind")
              + offCentre*fvc::div(phiBig, T, "implicit")
            );
            TEqn.solve();
        }
        divPhiOld = fvc::div(phiSmall, T, "explicit")
                  + fvc::div(phiBig, T, "implicit");
        
        scalarField TV = 0.25*fvc::surfaceSum
        (
            mag(fvc::snGrad(T))/mesh.deltaCoeffs()
        )().primitiveField();
        const scalarField& V = mesh.V().field();

        os << runTime.timeName() << " " << min(T.internalField()).value() << " "
           << max(T.internalField()).value() << " "
           << gSum(T.primitiveField()*V)/gSum(V) << " "
           << gSum(TV) << endl;

        runTime.write();

    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
