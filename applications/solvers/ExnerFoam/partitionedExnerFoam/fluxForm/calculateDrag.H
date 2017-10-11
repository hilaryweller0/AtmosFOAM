// Calculating the drag forced between two partitions

if (mag(dragCoeff) > SMALL)
{
    drag[0] = -fvc::interpolate(sigma[0]*sigma[1]*rhoSum)/cloudRadius*dragCoeff
            * mag(Uf[0] - Uf[1])*(Uf[0] - Uf[1]) & mesh.Sf();

    drag[1] = -drag[0];

    // The drag force shouldn't reverse the sign of the velocity. Check that 
    // drag*dt < flux
    Info << "Maximum dt*drag is " << dt.value()*max(mag(drag[0])).value() << nl
         << "Maximum stable.flux = " << max(mag(flux[0])) << nl
         << "Maximum buoyant.flux = " << max(mag(flux[1])) << endl;
}
