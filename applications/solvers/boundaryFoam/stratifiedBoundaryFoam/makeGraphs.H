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
volSymmTensorField S = dev(symm(gradU));
volScalarField S2(2*magSqr(S));
volScalarField magS(sqrt(S2));
volScalarField W
(
    (2*Foam::sqrt(2.0))*((S&S)&&S)
   /(
        magS*S2
      + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
    )
);
volScalarField phis
(
    (1.0/3.0)*acos(min(max(Foam::sqrt(6.0)*W, -scalar(1)), scalar(1)))
);
volScalarField As(Foam::sqrt(6.0)*cos(phis));
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
    "mixingLength",
    sqrt(turbulence->nut()/max(magS, dimensionedScalar("", dimless/dimTime, SMALL)))
);

volScalarField Richardson
(
    "Ri",
    dbdz/max(magSqr(vector(0,0,1) & gradU), 
            dimensionedScalar("", dbdz.dimensions(), VSMALL)
            )
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
makeGraph(z, W, "W", gFormat);
makeGraph(z, phis, "phi", gFormat);
makeGraph(z, As, "As", gFormat);
makeGraph(z, Us, "Us", gFormat);
makeGraph(z, C1, "C1", gFormat);
makeGraph(z, Richardson, "Ri", gFormat);

