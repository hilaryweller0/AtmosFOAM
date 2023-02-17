volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->R()
);

runTime.write();

// More diagnositcs to debug stratifiedRealizableKE
volTensorField gradU(fvc::grad(U));

surfaceTensorField gradUf = linearInterpolate(gradU);
gradUf += (fvc::snGrad(U)*mesh.magSf() - (mesh.Sf() & gradUf))
                *mesh.Sf()/sqr(mesh.magSf());
surfaceScalarField S2f(2*magSqr(dev(symm(gradUf))));
volScalarField S2 = fvc::average(S2f);

//volScalarField S2(2*magSqr(dev(symm(gradU))));
bound(S2, dimensionedScalar("", S2.dimensions(), VSMALL));
volScalarField magS(sqrt(S2));
volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

volScalarField eta
(
    "eta",
    magS*turbulence->k()/turbulence->epsilon()
);
volScalarField C1(max(eta/(scalar(5) + eta), scalar(0.43)));

// Diagnostics to write out

volScalarField mixingLength
(
    "mixingLength", sign(turbulence->nut())*sqrt(mag(turbulence->nut())/magS)
);

const word& gFormat = runTime.graphFormat();

makeGraph(z, U.component(0), "Ux", gFormat);
makeGraph(z, U.component(1), "Uy", gFormat);
makeGraph(z, turbulence->nut()/scalar(1), "nut", gFormat);
makeGraph(z, turbulence->k()/scalar(1), "k", gFormat);
makeGraph(z, turbulence->epsilon()/scalar(1), "epsilon", gFormat);
makeGraph(z, magS, "magS", gFormat);
makeGraph(z, eta, gFormat);
makeGraph(z, mixingLength, "lm", gFormat);
makeGraph(z, mesh.lookupObject<const volScalarField>("Cmu"), "Cmu", gFormat);
makeGraph(z, mag(gradU), "magGradU", gFormat);
makeGraph(z, Us, "Us", gFormat);
makeGraph(z, C1, "C1", gFormat);
makeGraph(z, mesh.lookupObject<const volScalarField>("Ri"), "Ri", gFormat);

