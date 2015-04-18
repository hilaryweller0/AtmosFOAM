/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Application
    meshAnalysis2D

Description
    Calculates the non-orthogonality, skewness and some other grid
    characteristics and writes the results relative to distance
    from a given point

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellQuality.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::validArgs.append("dict");
    #include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    runTime.setTime(Times[startTime], startTime);

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    const fileName dictName(args.additionalArgs()[0]);
    
    // Read the centre point for the analysis from the dictionary
    IOdictionary dict
    (
        IOobject(dictName, runTime.system(), runTime, IOobject::MUST_READ)
    );
    const point centre(dict.lookup("centre"));

    cellQuality cq(mesh);

    surfaceScalarField orthogonality
    (
        IOobject("orthogonality", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("degrees", dimless, scalar(0))
    );
    const scalarField orthogTmp = cq.faceNonOrthogonality();
    forAll(orthogonality, faceI)
    {
        orthogonality[faceI] = orthogTmp[faceI];
    }
    orthogonality.write();
    
    surfaceScalarField skewness
    (
        IOobject("skewness", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("degrees", dimless, scalar(0))
    );
    const scalarField skewTmp = cq.faceSkewness();
    forAll(skewness, faceI)
    {
        skewness[faceI] = skewTmp[faceI];
    }
    skewness.write();
        
    surfaceScalarField dx("dx", 1./mesh.deltaCoeffs());
    dx.write();
    
    const dimensionedScalar Vtot = sum(mesh.V());

    Info << "n cells = " << mesh.nCells()
         << " ndofs = " << mesh.nCells()+mesh.nInternalFaces()
         << " max/min dx = " << (max(dx)/min(dx)).value()
         << " max dx = " << (max(dx)).value()
         << "\nmax non-orthogonality = " << max(orthogonality).value()
         << "\nmax skewness = " << max(skewness).value() << endl;
    
    // Read detHess
    volScalarField detHess
    (
        IOobject("detHess", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("detHess", dimless, scalar(0))
    );

    // Write the mesh analysis for volFields
    fileName ofName
    (
        runTime.rootPath()/runTime.caseName()/runTime.timeName()/"distArea.dat"
    );
    Info << "Writing to " << ofName << endl;
    OFstream ofVol(ofName);
    ofVol << "#dist area detHess\n";
    forAll(mesh.V(), cellI)
    {
        scalar dist = mag(mesh.C()[cellI] - centre);
        ofVol << dist << " " << mesh.V()[cellI] << " " << detHess[cellI] << nl;
    }
    
    // Write the mesh analysis for surfaceFields
    fileName ofName2
    (
        runTime.rootPath()/runTime.caseName()/runTime.timeName()/"distMetrics.dat"
    );
    Info << "Writing to " << ofName2 << endl;
    OFstream ofS(ofName2);
    ofS << "#dist dx orthogonality skewness\n";
    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        scalar dist = mag(mesh.Cf()[faceI] - centre);
        ofS << dist << " " << dx[faceI] << " " << orthogonality[faceI]
            << " " << skewness[faceI] << nl;
    }
    
    // Calculat the scalaring factors between the monitor fucntion and the cell areas
    // if monitor function given
    if (dict.found("monitorFunctionFrom"))
    {
        // Budd et al monitor funciton
        if (word(dict.lookup("monitorFunctionFrom")) == "monitorFunctionSech")
        {
            scalar alpha1 = readScalar(dict.lookup("alpha1"));
            scalar alpha2 = readScalar(dict.lookup("alpha2"));
            scalar a  = readScalar(dict.lookup("a"));
        
            scalar sum = 0;
            forAll(mesh.C(), cellI)
            {
                scalar dist = mag(mesh.C()[cellI] - centre);
                scalar rho = 1 + alpha1/sqr(Foam::cosh(alpha2*(sqr(dist) - sqr(a))));
                sum += 1/rho;
            }
            
            const scalar theoryScale = 1./sum;
            Info << "Scaling factor between theoretical cell area and achieved cell "
                 << "area = " << theoryScale << endl;
            
            // Output the theoretical dx and cell area
            OFstream ofVol
            (
                runTime.rootPath()/runTime.caseName()
               /runTime.timeName()/"distAreaTheory.dat"
            );
            ofVol << "#dist dx area" << endl;
            forAll(mesh.C(), cellI)
            {
                scalar dist = mag(mesh.C()[cellI] - centre);
                scalar rho = 1 + alpha1/sqr(Foam::cosh(alpha2*(sqr(dist) - sqr(a))));
                scalar area = theoryScale/rho;
                // Cell-centre to cell centre distance for a square
                scalar dx = Foam::sqrt(area);
                ofVol << dist << " " << dx << " " << area << nl;
            }
        }
       // Tanh monitor funciton from Ringler et al
        else if (word(dict.lookup("monitorFunctionFrom")) == "monitorFunctionTanh")
        {
            scalar alpha = readScalar(dict.lookup("alpha"));
            scalar beta = readScalar(dict.lookup("beta"));
            scalar gamma = readScalar(dict.lookup("gamma"));
    
            scalar sum = 0;
            forAll(mesh.C(), cellI)
            {
                scalar dist = mag(mesh.C()[cellI] - centre);
                scalar rho = 0.5/(1+gamma)*(Foam::tanh((beta-dist)/alpha)+1) + gamma;
                sum += 1/Foam::sqrt(rho);
            }
            
            const scalar theoryScale = 1./sum;
            Info << "Scaling factor between theoretical cell area and achieved cell "
                 << "area = " << theoryScale << endl;
            
            // Output the theoretical dx and cell area
            OFstream ofVol
            (
                runTime.rootPath()/runTime.caseName()
               /runTime.timeName()/"distAreaTheory.dat"
            );
            ofVol << "#dist dx area" << endl;
            forAll(mesh.C(), cellI)
            {
                scalar dist = mag(mesh.C()[cellI] - centre);
                scalar rho = 0.5/(1+gamma)*(Foam::tanh((beta-dist)/alpha)+1) + gamma;
                scalar area = theoryScale/Foam::sqrt(rho);
                // Cell-centre to cell centre distance for a square
                scalar dx = Foam::sqrt(area);
                ofVol << dist << " " << dx << " " << area << nl;
            }
        }
    }
}


// ************************************************************************* //


