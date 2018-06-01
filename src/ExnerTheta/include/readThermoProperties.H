Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar T0
(
    "T0", dimTemperature, readScalar(thermoDict.lookup("T0"))
);
const dimensionedScalar pRef
(
    "pRef", dimPressure, readScalar(thermoDict.lookup("pRef"))
);

const constTransport<hConstThermo<perfectGas<specie> > > air
(
    thermoDict.subDict("mixture")
);

const dimensionedScalar R("R", dimGasConstant, air.R());
const dimensionedScalar Cp
(
    "Cp", dimGasConstant, air.Cp(pRef.value(),T0.value())
);
const dimensionedScalar kappa("kappa", R/Cp);
const scalar gamma = (1/(1-kappa)).value();
const dimensionedScalar Cv("Cv", Cp - R);
const dimensionedScalar speedSound("speedSound", sqrt(gamma*R*T0));

Info << "Cp = " << Cp.value() << " speedSound = " << speedSound << " kappa = "
     << kappa << " R = " << air.R() << " RR = " << RR << " Cv = " << Cv << endl;

Info << "Maximum acousitic Courant number = "
     << max(mesh.deltaCoeffs())*runTime.deltaT()*speedSound << endl;

