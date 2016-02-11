/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    createGhostMesh

Description
    Takes a mesh with cyclic boundary conditions and extends it, creating ghost
    cells adjacent to the cyclic boundaries which are the same as the cells on
    adjacent to the cyclic neighbour patch. This can then be used for
    mapToGhost to map values from the original to the ghost mesh, including 
    ghost cells. Then you can run a tailored application. Cyclic boundaries
    are treated in this way in order to be able to use the extended stencil
    differencing schemes
    
    The argument is the number of overlap cells in each ghost region
    
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "Time.H"
#include "argList.H"
#include "emptyFvPatch.H"
#include "cyclicFvPatch.H"
#include "fvMeshSubset.H"
#include "mergePolyMesh.H"
#include "polyTopoChanger.H"
#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Create a mesh with ghost cells instead of cyclic BCs");
    argList::validArgs.append("nOverlap");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The number of overlap cells
    label nOverlap = readLabel(IStringStream(args[1])());
    
    Info << "Original mesh has " << mesh.nCells()
         << " cells. Adding overlap regions with " << nOverlap
         << " cells to each cyclic boundary" << endl;

    // Reference to the patches of the original mesh
    const fvPatchList& oldFvPatches = mesh.boundary();
    const polyPatchList& oldPatches = mesh.boundaryMesh();

    // List of cyclic patches and their neighbours
    labelList cyclics(oldFvPatches.size());
    labelList cyclicNeighbs(oldFvPatches.size());
    wordList cyclicNames(oldFvPatches.size());
    wordList cyclicNeighbours(oldFvPatches.size());
    wordList cyclicSubMeshNames(oldFvPatches.size());
    List<vector> cyclicDeltas(oldFvPatches.size());
    label nCyclics = 0;

    // First find the cyclic boundaries and the transforms between cyclics
    for(label ip = 0; ip < oldFvPatches.size(); ip++)
    {
        if (isA<cyclicFvPatch>(oldFvPatches[ip]))
        {
            cyclics[nCyclics] = ip;
        
            // Reference to the cyclicFvPatch
            const cyclicFvPatch& oldFvPatch
                 = refCast<const cyclicFvPatch>(oldFvPatches[ip]);
        
            // Store the names of the cyclic BCs and their neighbour Patches
            cyclicNames[nCyclics] = oldFvPatch.name();
            cyclicNeighbs[nCyclics] = oldFvPatch.neighbPatchID();
            
            // Only works for parallel cyclics with one transformation vector
            if (!oldFvPatch.parallel())
            {
                FatalErrorIn("createGhostMesh")
                    << "only works for parallel cyclic boundaries but "
                    << oldFvPatch.name() << " and its neighbour are not parallel"
                    << exit(FatalError);
            }

            // Find the constant vector across the cyclic pair
            const pointField& patchPoints = oldFvPatch.patch().localPoints();
            pointField transPoints(1, patchPoints[0]);
            oldFvPatch.cyclicPatch().transformPosition(transPoints);
            cyclicDeltas[nCyclics] = patchPoints[0] - transPoints[0];

            nCyclics++;
        }
    }
    
    // Replace the cyclic boundaries with patches in the original mesh
    List<polyPatch*> newPatches(mesh.boundaryMesh().size());

    for(label ip = 0; ip < oldPatches.size(); ip++)
    {
        if (isA<emptyPolyPatch>(oldPatches[ip]))
        {
            newPatches[ip] = new emptyPolyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                oldPatches[ip].start(), ip, mesh.boundaryMesh(), "empty"
            );
        }
        else
        {
            newPatches[ip] = new polyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                oldPatches[ip].start(), ip, mesh.boundaryMesh(), "patch"
            );
        }
    }
    mesh.removeFvBoundary();
    mesh.addFvPatches(newPatches);

    // Create the ghostMesh as a mergePolyMesh so that subMeshes can be added
    mergePolyMesh ghostMesh
    (
        IOobject
        (
            fvMesh::defaultRegion, runTime.timeName(), runTime,
            IOobject::MUST_READ
        )
    );
    ghostMesh.rename("ghostMesh");
    ghostMesh.setInstance("ghostMesh");
    
    // Maps from mesh to ghostMesh cells
    IOList<label> meshToGhostCells
    (
        IOobject
        (
            "ghostToMeshCells",
            ghostMesh.time().timeName(),
            ghostMesh.meshSubDir,
            ghostMesh
        ),
        labelList(0)
    );
    
    // Add boundaries to the ghostMesh
    List<polyPatch*> ghostPatches(newPatches.size());
    for(label ip = 0; ip < newPatches.size(); ip++)
    {
        if (isA<emptyPolyPatch>(oldPatches[ip]))
        {
            ghostPatches[ip] = new emptyPolyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                oldPatches[ip].start(), ip, ghostMesh.boundaryMesh(), "empty"
            );
        }
        else
        {
            ghostPatches[ip] = new polyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                oldPatches[ip].start(), ip, ghostMesh.boundaryMesh(), "patch"
            );
        }
    }
    ghostMesh.removeBoundary();
    ghostMesh.addPatches(ghostPatches);
    
    // Add sub meshes to ghostMesh for each cyclic boundary
    for(label ipc = 0; ipc < nCyclics; ipc++)
    {
        label ip = cyclics[ipc];
        
        const string ipStr = static_cast<std::ostringstream*>
                ( &(std::ostringstream() << ip) )->str();
        const string ipStrn = static_cast<std::ostringstream*>
                ( &(std::ostringstream() << cyclicNeighbs[ipc]) )->str();
        cyclicSubMeshNames[ipc]
             = mesh.boundary()[cyclicNeighbs[ipc]].name()+ipStrn;

        const fvPatch& oldFvPatch = mesh.boundary()[ip];
            
        // Find the cells adjacent to the cyclic boundary
        labelHashSet bCells;
        
        if (nOverlap >= 1)
        {
            const labelList& faceCells = oldFvPatch.faceCells();
            bCells.insert(faceCells);
        }

        // Add more cells to the submesh if necessary
        for(label iover = 1; iover < nOverlap; iover++)
        {
            labelList bCellsOld(bCells.sortedToc());
            forAll(bCellsOld, cellI)
            {
                bCells.insert(mesh.cellCells()[bCellsOld[cellI]]);
            }
        }

        // Create the submesh based on the cells adjacent to the boundary
        fvMeshSubset subsetMesh(mesh);
        subsetMesh.setCellSubset(bCells, cyclicNeighbs[ipc]);
            
        // Manipulate the sub-mesh
        fvMesh& subMesh = subsetMesh.subMesh();
            
        // Moves the points of the subMesh by the cyclic transformation
        pointField newPoints(subMesh.points()+cyclicDeltas[ipc]);
        subMesh.movePoints(newPoints);
            
        // Rename patches of the submesh
        const polyPatchList& oldSubPatches = subMesh.boundaryMesh();
        List<polyPatch*> subPatches(oldSubPatches.size());
        for(label ipp = 0; ipp < oldSubPatches.size(); ipp++)
        {
            if (isA<emptyPolyPatch>(oldSubPatches[ipp]))
            {
                subPatches[ipp] = new emptyPolyPatch
                (
                    oldSubPatches[ipp].name()+ipStr, oldSubPatches[ipp].size(),
                    oldSubPatches[ipp].start(), ipp, subMesh.boundaryMesh(),
                    "empty"
                );
            }
            else
            {
                subPatches[ipp] = new polyPatch
                (
                    oldSubPatches[ipp].name()+ipStr, oldSubPatches[ipp].size(),
                    oldSubPatches[ipp].start(), ipp, subMesh.boundaryMesh(),
                    "patch"
                );
            }
        }
        subMesh.removeFvBoundary();
        subMesh.addPatches(subPatches);

        // Add the subMesh to the ghostMesh
        Info << "Adding " << subMesh.nCells() << endl;
        ghostMesh.addMesh(subMesh);

        // Add the cell numbers to the map for mapping to ghost cells
        meshToGhostCells.append(bCells.sortedToc());
    }

    ghostMesh.merge();
    Info << "Writing final ghostMesh with " << ghostMesh.nCells()
         << " cells" << endl;
    ghostMesh.write();
    meshToGhostCells.write();
    
    Info << "Run:\n";
    for(label ip = 0; ip < nCyclics; ip++)
    {
        Info << "stitchMesh -perfect -overwrite -region ghostMesh "
             << cyclicNames[ip] << " " << cyclicSubMeshNames[ip] << endl;
    }
}
