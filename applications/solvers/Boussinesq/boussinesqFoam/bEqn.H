// Advective form buoyancy equation
dbdt = - fvc::div(volFlux, bBar + b, "b") 
              + (bBar + b)*fvc::div(volFlux);

// Viscous term
if (alpha.value() > VSMALL)
{
    dbdt += alpha * fvc::laplacian(b);
}

if (!adjustQ)
{
    dbdt += Q;
}

if (relaxbTime.value() > 0)
{
    dbdt -= b/relaxbTime;
}

// Update b according to the rate of change
b = b.oldTime() + dt*
(
    (1-offCentre)*dbdt.oldTime()
  + offCentre*dbdt
);
b.correctBoundaryConditions();

// bouyancy force on faces
bf = fvc::interpolate(b,"b")*mesh.Sf().component(2);

if (adjustQ)
{
    dimensionedScalar totalHeat("totalHeat", dimensionSet(0, 1,-2,0,0), scalar(0));
    totalHeat = fvc::domainIntegrate(b)/sum(mesh.V());
    Q = (totalHeat0 - totalHeat)/dt;
    b += dt*Q;
}

