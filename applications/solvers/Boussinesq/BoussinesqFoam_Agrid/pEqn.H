// Pressure equation
U = UEqn.H()/UEqn.A();
surfaceScalarField phiTmp = fvc::interpolate(U) & mesh.Sf();
gradPcoeff = fvc::interpolate(1/UEqn.A());

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(gradPcoeff, p) ==  fvc::div(phiTmp) + fvc::div(bf)
    );
    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();
    
    if(nonOrth == nNonOrthCorr)
    {
        phi = phiTmp + bf - pEqn.flux();
    }
}

U += (b*kdir - fvc::grad(p))/UEqn.A();
U.correctBoundaryConditions();
