// Calculating the drag forced between two partitions

if (mag(dragCoeff) > SMALL)
{
    PtrList<surfaceScalarField> sigmaf(nParts);
    for (label ip = 0; ip < nParts; ip++)
    {
        sigmaf.set
        (
            ip,
            new surfaceScalarField(fvc::interpolate(sigma[ip]))
        );
    }
    
    surfaceScalarField cloudRadius
    (
        "cloudRadius",
        max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1])
    );
    surfaceScalarField dragCommon = -dragCoeff*fvc::interpolate(rhoSum)
                 /cloudRadius * mag(Uf[0] - Uf[1])*(Uf[0] - Uf[1]) & mesh.Sf();

    drag[0] = sigmaf[1]/fvc::interpolate(rho[0])*dragCommon;
    drag[1] = -sigmaf[0]/fvc::interpolate(rho[1])*dragCommon;

    surfaceScalarField netDrag = drag[0]*sigmaf[0]*fvc::interpolate(rho[0])
                               + drag[1]*sigmaf[1]*fvc::interpolate(rho[1]);

    // The drag force shouldn't reverse the sign of the velocity. Check that 
    // drag*dt < flux
    Info << "Maximum dt*drag is " << dt.value()*max(mag(drag[0])).value() << nl
         << "Maximum stable.flux = " << max(mag(volFlux[0])) << nl
         << "Maximum buoyant.flux = " << max(mag(volFlux[1])) << endl;
}
