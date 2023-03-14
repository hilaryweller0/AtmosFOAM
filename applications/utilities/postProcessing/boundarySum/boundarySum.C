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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    boundarySum

Description
    Calculates various sums over boundary patches.
    Output is in a file called "sum_"+boundaryPatchName+fieldName+".dat. 
    The  sums are for the input variable T = fieldName
    Each line
    of the output file consists of 
    time  L1 L2 Linf L0 variance min max
    where L1 = 1/S integral(|T| dS)
          L2 = 1/S integral(T^2 dS)
          Li = max(|T|) over the surface
          L0 = 1/S int(T dS)
          variance = 1/S int((T - L0)^2 dS)
          min = min(T) over the domain
          max = max(T) over the domain
    where S is the area of the boundary patch

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("field");
    argList::validArgs.append("boundaryName");

    argList::addNote
    (
        "Calculates various sums over a boundary patch.\n"
        "Output is in a file called sum_+boundaryPatchName+fieldName+.dat or\n"
        "The sums are for the input variable T = fieldName\n"
        "Each line\n"
        "of the output file consists of\n"
        "time  L1 L2 Linf L0 variance min max\n"
        "where L1 = 1/S integral(|T| dS)\n"
        "      L2 = 1/S integral(T^2 dS)\n"
        "      Li = max(|T|) over the boundary patch\n"
        "      L0 = 1/S int(T dS)\n"
        "      variance = 1/S int((T - L0)^2 dS)\n"
        "      min = min(T) over the patch\n"
        "      max = max(T) over the patch\n"
        "where S is the area of the boundary patch\n"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    const word fieldName(args.args()[1].c_str());
    const word boundaryName(args.args()[2].c_str());

    #include "createMesh.H"
    
    // Find the boundary patch
    const label ipatch = mesh.boundary().findPatchID(boundaryName);

    // initialise diagnostics file
    fileName outFile = 
        args.rootPath() / args.caseName() /"sum_"+boundaryName+fieldName+".dat";
    Info << "Writing boundary sums to " << outFile << endl;
    OFstream os(outFile);
    os << "#time mag RMS inf sum variance min max" << endl;
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        mesh.readUpdate();

        IOobject header
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (!header.headerOk())
        {
            Info << "Time = " << runTime.timeName() << " no " << fieldName
                 << endl;
        }
        else if(header.headerClassName() == "volScalarField")
        {
            const volScalarField vf(header, mesh);
            const fvPatchField<scalar>& f = vf.boundaryField()[ipatch];
            
            scalar Stot = 0;
            scalar l1 = 0;
            scalar l2 = 0;
            scalar li = 0;
            scalar l0 = 0;
            scalar min = GREAT;
            scalar max = -GREAT;
            scalar variance = 0;
            
            forAll(f, i)
            {
                scalar fi = f[i];
                scalar Si = mesh.boundary()[ipatch].magSf()[i];
                Stot += Si;
                l1 += mag(fi)*Si;
                l2 += sqr(fi)*Si;
                if (mag(fi) > li) li = mag(fi);
                l0 += fi*Si;
                variance += sqr(fi)*Si;
                if (fi > max) max = fi;
                if (fi < min) min = fi;
            }
            l1 /= Stot;
            l2 = Foam::sqrt(l2/Stot);
            l0 /= Stot;
            variance = (variance/Stot - sqr(l0));

            os << runTime.timeName() << ' ' << l1 << ' ' << l2 << ' ' << li
               << ' ' << l0 << ' ' << variance << ' ' << min << ' ' << max << endl;
        }
        else
        {
            FatalErrorIn("boundarySum") << "Field type "
                << header.headerClassName() << " of " << fieldName
                << " not supported" << exit(FatalError);
        }
    }
    return(0);
}


// ************************************************************************* //
