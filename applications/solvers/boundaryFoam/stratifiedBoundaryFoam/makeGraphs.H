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
    turbulence->sigma()
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

setWriter::New(runTime.graphFormat())->write
(
    runTime.globalPath()
   /functionObjects::writeFile::outputPrefix
   /args.executable()
   /runTime.timeName(),

    args.executable(),

    coordSet(true, "z", z),

    "Ux", (U.component(0))().primitiveField(),
    "Uy", (U.component(1))().primitiveField(),

    "nu", turbulence->nu()().primitiveField(),
    "k", turbulence->nu()().primitiveField(),
    "epsilon", turbulence->nu()().primitiveField(),

    "magS", magS.primitiveField(),
    "eta", eta.primitiveField(),
    "Cmu", (mesh.lookupObject<const volScalarField>("Cmu")).primitiveField(),

    "gammaDot", mag(fvc::grad(U))().primitiveField(),
    "lm", mixingLength.primitiveField(),

    "Us", Us.primitiveField(),
    "C1", C1.primitiveField(),
    "Ri", mesh.lookupObject<const volScalarField>("Ri").primitiveField(),
    "dbdz", dbdz.primitiveField()
);

