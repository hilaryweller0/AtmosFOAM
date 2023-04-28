{
    // Reading in the sponge layer coefficients (if present)
    surfaceScalarField muSponge
    (
        IOobject("muSponge", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("muSponge", dimless, scalar(0))
    );
    Info << "Maximum sponge value = " << max(muSponge).value() << endl;

    vectorField& faceAreas = const_cast<vectorField&>(mesh.faceAreas());
    scalarField& magFaceAreas = const_cast<scalarField&>(mesh.magFaceAreas());

    for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        faceAreas[facei]    *= 1-muSponge[facei];
        magFaceAreas[facei] *= 1-muSponge[facei];
    }
}
