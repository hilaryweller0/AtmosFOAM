{
    // Evaluate bulk behaviour
    scalar delta = mesh.bounds().span()[1];
    scalar nu = turbulence->nu()().boundaryField()[patchId][faceId];
    scalar Rey = Ubar[0].value()*delta/nu;
    symmTensor R = turbulence->devReff()().boundaryField()[patchId][faceId];
    scalar tauw = flowDirection & R & wallNormal;
    scalar uTau = ::sqrt(mag(tauw));
    scalar ReyFric = uTau*delta/nu;

    Info<< "Rey = " << Rey << ", ReyFric = " << ReyFric
        << endl;
}
