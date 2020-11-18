fvScalarMatrix bEqn
(
    fvm::ddt(b)
  + fvm::div(phi, b)
  - fvm::laplacian(alpha, b)
);

if (implicitAdvection)
{
    solve(bEqn == Q);
}
else
{
    b = (bEqn.H() + Q)/bEqn.A();
    b.correctBoundaryConditions();
    b = (bEqn.H() + Q)/bEqn.A();
    b.correctBoundaryConditions();
}

// Bouyancy on the faces in the z direction
bf = fvc::interpolate(b)*mesh.Sf().component(2)/fvc::interpolate(UEqn.A());
