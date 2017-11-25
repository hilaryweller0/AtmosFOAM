// Calculating the drag forced between two partitions

if (mag(dragCoeff) > SMALL)
{
    surfaceScalarField cloudRadius
    (
        "cloudRadius",
        min
        (
            fvc::interpolate(sigma[0])*fvc::interpolate(sigma[1])*cloudRadiusMin,
            cloudRadiusMax
        )
    );

    drag[0] = -dragCoeff/cloudRadius * mag(Uf[0] - Uf[1])*(Uf[0] - Uf[1]) & mesh.Sf();

    drag[0] = fvc::interpolate(sigma[1])*drag[0];
    drag[1] = -fvc::interpolate(sigma[0])*drag[0]/fvc::interpolate(sigma[1])
              *fvc::interpolate(rho[0])/fvc::interpolate(rho[1]);

    // The drag force shouldn't reverse the sign of the velocity. Check that 
    // drag*dt < flux
    Info << "Maximum dt*drag is " << dt.value()*max(mag(drag[0])).value() << nl
         << "Maximum stable.flux = " << max(mag(volFlux[0])) << nl
         << "Maximum buoyant.flux = " << max(mag(volFlux[1])) << endl;
}
