if (divFreeInitial)
{
    Info << "Making the initial conditions divergence free" << endl;

    volScalarField P
    (
        IOobject("P", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimensionSet(1, -1,-1,0,0), scalar(0)),
        "zeroGradient"
    );

    fvScalarMatrix PEqn
    (
        fvc::div(phi) - fvm::laplacian(P)
    );
    PEqn.setReference(0,0);
    PEqn.solve("ExnerFinal");
    phi += PEqn.flux();
}
