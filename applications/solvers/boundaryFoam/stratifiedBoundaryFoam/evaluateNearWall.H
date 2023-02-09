{
    // Evaluate near-wall behaviour

    scalar nu = turbulence->nu()().boundaryField()[patchId][faceId];
    scalar nut = turbulence->nut()().boundaryField()[patchId][faceId];
    symmTensor R = turbulence->devReff()().boundaryField()[patchId][faceId];
    scalar epsilon = turbulence->epsilon()()[cellId];
//    scalar omega = turbulence->omega()()[cellId];
    scalar k = turbulence->k()()[cellId];
    scalar magUp = mag(U[cellId] - U.boundaryField()[patchId][faceId]);

    vector flowDirection = U[cellId]/mag(U[cellId]);
    scalar tauw = flowDirection & R & wallNormal;

    scalar uTau = ::sqrt(mag(tauw));

    scalar zPlus = uTau*z[cellId]/(nu + rootVSmall);

    scalar uPlus = magUp/(uTau + rootVSmall);

    scalar nutPlus = nut/nu;

    scalar kPlus = k/(sqr(uTau) + rootVSmall);

    scalar epsilonPlus = epsilon*nu/(pow4(uTau) + rootVSmall);

//    scalar omegaPlus = omega*nu/(sqr(uTau) + rootVSmall);

    scalar Rey = magUp*z[cellId]/nu;

    Info<< runTime.timeName() << " Rey = " << Rey << " uTau = " << uTau << " nut+ = " << nutPlus
        << " z+ = " << zPlus << " u+ = " << uPlus
        << " k+ = " << kPlus << " epsilon+ = " << epsilonPlus
        << endl;
}
