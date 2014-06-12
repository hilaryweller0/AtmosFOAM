//
// Make the boundaries of the fvMesh orthogonal (by moving boundary face centres)
// ~~~~~~~~~~~~

    Foam::Info << "Making the boundary faces orthogonal" << Foam::endl;

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& pat = mesh.boundary()[patchi];
        // only move centres of un-coupled boundary faces
        if (!pat.coupled())
        {
            // loop through and move all of the faces
            for(label facei = pat.start(); facei <pat.start()+pat.size();facei++)
            {
                label own = mesh.faceOwner()[facei];
                const vector& C = mesh.C()[own];
                const vector& Cf = mesh.faceCentres()[facei];
                const vector Sfhat = mesh.faceAreas()[facei]
                                      /mag(mesh.faceAreas()[facei]);
                const vector fCnew = C + ((Cf - C) & Sfhat)*Sfhat;

                const_cast<vectorField&>(mesh.faceCentres())[facei] = fCnew;
            }
        }
    }
