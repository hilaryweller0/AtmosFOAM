// Ensure that empty faces are in the final patch of the ghostMesh
{
    const polyPatchList& oldPatches = ghostMesh.boundaryMesh();
    List<polyPatch*> newPatches(oldPatches.size());
    
    label ipp = 0;
    label newStart = oldPatches[0].start();
    // First add non-empty patches
    for(label ip = 0; ip < oldPatches.size(); ip++)
    {
        if (!isA<emptyPolyPatch>(oldPatches[ip]))
        {
            ghostPatches[ip] = new polyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                newStart, ipp, ghostMesh.boundaryMesh(), "patch"
            );
            newStart += oldPatches[ip].size();
            ipp++;
        }
    }

    // Next add empty patches
    for(label ip = 0; ip < oldPatches.size(); ip++)
    {
        if (isA<emptyPolyPatch>(oldPatches[ip]))
        {
            ghostPatches[ip] = new emptyPolyPatch
            (
                oldPatches[ip].name(), oldPatches[ip].size(),
                newStart, ipp, ghostMesh.boundaryMesh(), "empty"
            );
            newStart += oldPatches[ip].size();
            ipp++;
        }
    }

    ghostMesh.removeBoundary();
    ghostMesh.addPatches(newPatches);
}
