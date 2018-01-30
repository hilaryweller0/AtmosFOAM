// Calculating the drag force between two partitions

if (mag(dragCoeff) > SMALL)
{
    // This code only works for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("calculateDrag.H")
            << " drag is only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    // Calculate the cloud radius and hence the drag
    cloudRadius = max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1]);
    dragCommon = -dragCoeff/cloudRadius 
                 * mag(Uf[0] - Uf[1])*(Uf[0] - Uf[1]) & mesh.Sf();

    drag[0] = sigmaf[1]*dragCommon*sigmaRhof.sum()
                /fvc::interpolate(rho[0], "rho");
    drag[1] =-sigmaf[0]*dragCommon*sigmaRhof.sum()
                /fvc::interpolate(rho[1], "rho");

    // The drag force shouldn't reverse the sign of the velocity. Check that 
    // drag*dt < flux
    Info << "Maximum dt*drag[0] is " << dt.value()*max(mag(drag[0])).value() << nl
         << "Maximum dt*drag[1] is " << dt.value()*max(mag(drag[1])).value() << nl
         << "Maximum stable.flux = " << max(mag(volFlux[0])) << nl
         << "Maximum buoyant.flux = " << max(mag(volFlux[1])) << endl;
}
